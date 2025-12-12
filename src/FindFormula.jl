module FindFormulaModule

# Access ELEMENTS from parent module
const ELEMENTS = if isdefined(parentmodule(@__MODULE__), :ELEMENTS)
    parentmodule(@__MODULE__).ELEMENTS
else
    error("ELEMENTS constant not found in parent module. Ensure IsotopicCalc is properly initialized with elements.json data.")
end

"""
    Compound

Data structure representing a molecular formula match from mass spectrometry data.

# Fields
- `formula::String`: Molecular formula in Hill notation (e.g., "C3H6O")
- `adduct::String`: Adduct display name (e.g., "M+H", "M+Na", or "" for none)
- `charge::Int`: Charge state (positive, negative, or 0 for neutral)
- `mz::Float64`: Calculated m/z value in atomic mass units
- `ppm::Float64`: Mass accuracy in parts per million (ppm)

# Examples
```julia
julia> compound = Compound("C3H6O", "M+H", 1, 59.049141, -0.7)
Compound("C3H6O", "M+H", 1, 59.049141, -0.7)
```

# Notes
The `adduct` field stores the display name (e.g., "M+H") for output formatting,
while the actual adduct information (mass, charge) is handled by the `AdductInfo` struct.

See also: [`find_formula`](@ref), [`AdductInfo`](@ref)
"""
struct Compound
    formula::String
    adduct::String
    charge::Int
    mz::Float64
    ppm::Float64
end

"""
    AdductInfo

Data structure representing an adduct type with its mass and charge.

# Fields
- `mass::Float64`: Mass change in atomic mass units (amu)
- `charge::Int`: Charge state (positive, negative, or 0 for neutral)
- `name::String`: Display name for the adduct
"""
struct AdductInfo
    mass::Float64
    charge::Int
    name::String
end

"""
Dictionary of common mass spectrometry adducts.

# Supported Adducts
- "H+": Protonated [M+H]+ (+1.00783 amu, charge +1)
- "Na+": Sodium adduct [M+Na]+ (+22.98977 amu, charge +1)
- "K+": Potassium adduct [M+K]+ (+38.96371 amu, charge +1)
- "H-": Deprotonated [M-H]- (-1.00783 amu, charge -1)
- "+": Radical cation [M]+• (charge +1, electron mass correction applied during calculation)
- "": Neutral molecule (0 amu, charge 0)

# Notes
For charged species, electron mass correction (-0.0005485 amu per positive charge,
+0.0005485 amu per negative charge) is automatically applied during m/z calculation.
The radical cation "+" has no additional mass shift beyond the electron correction.
"""
const ADDUCTS = Dict(
    "H+" => AdductInfo(1.00782503223, 1, "M+H"),
    "Na+" => AdductInfo(22.989769282, 1, "M+Na"),
    "K+" => AdductInfo(38.9637064864, 1, "M+K"),
    "H-" => AdductInfo(-1.00782503223, -1, "M-H"),
    "+" => AdductInfo(0.0, 1, "M"),
    "" => AdductInfo(0.0, 0, "")
)

"""
    get_adduct_info(adduct::String) -> AdductInfo

Get adduct information from an adduct string.

# Arguments
- `adduct::String`: Adduct identifier (e.g., "H+", "Na+", "K+", "H-", or "")

# Returns
- `AdductInfo`: Struct containing mass, charge, and display name

# Throws
- `ArgumentError`: If the adduct string is not recognized

# Examples
```julia
julia> get_adduct_info("H+")
AdductInfo(1.00782503223, 1, "M+H")

julia> get_adduct_info("Na+")
AdductInfo(22.989769282, 1, "M+Na")
```
"""
function get_adduct_info(adduct::String)
    if !haskey(ADDUCTS, adduct)
        valid_adducts = join(sort(["\"$k\"" for k in keys(ADDUCTS)]), ", ")
        throw(ArgumentError("Unknown adduct '$adduct'. Supported adducts: $valid_adducts"))
    end
    return ADDUCTS[adduct]
end

"""
    build_hill_notation(elements, current) -> String

Construct a formula string in Hill notation order.

Hill notation rules:
- With carbon: C first, H second, then others alphabetically
- Without carbon: All elements alphabetically (including H)

# Arguments
- `elements::Vector{String}`: Element symbols
- `current::Vector{Int}`: Corresponding element counts

# Returns
- `String`: Formula in Hill notation (e.g., "C3H6O")
"""
function build_hill_notation(elements::Vector{String}, current::Vector{Int})
    hill_elements = String[]
    hill_counts = Int[]
    c_idx = findfirst(isequal("C"), elements)
    h_idx = findfirst(isequal("H"), elements)

    # Check if carbon is present
    has_carbon = !isnothing(c_idx) && current[c_idx] > 0

    if has_carbon
        # Carbon present: C, then H, then others alphabetically
        push!(hill_elements, "C")
        push!(hill_counts, current[c_idx])
        if !isnothing(h_idx) && current[h_idx] > 0
            push!(hill_elements, "H")
            push!(hill_counts, current[h_idx])
        end
        # Add remaining elements (excluding C and H) in alphabetical order
        remaining_idxs = [i for i in 1:length(elements) if current[i] > 0 && elements[i] != "C" && elements[i] != "H"]
        sorted_remaining = sort([(elements[i], current[i]) for i in remaining_idxs], by=x->x[1])
        for (el, cnt) in sorted_remaining
            push!(hill_elements, el)
            push!(hill_counts, cnt)
        end
    else
        # No carbon: All elements alphabetically (including H)
        all_idxs = [i for i in 1:length(elements) if current[i] > 0]
        sorted_all = sort([(elements[i], current[i]) for i in all_idxs], by=x->x[1])
        for (el, cnt) in sorted_all
            push!(hill_elements, el)
            push!(hill_counts, cnt)
        end
    end

    return join([string(hill_elements[i], hill_counts[i] == 1 ? "" : hill_counts[i]) for i in 1:length(hill_elements)])
end

"""
    generate_formulas(mz_input, atom_pool, adduct_info) -> Vector{Compound}

Generate all possible molecular formulas within the atom pool constraints.

Uses combinatorial search with chemical heuristics to enumerate formulas.
Currently applies H ≥ 0.5×C rule to reduce search space.

# Arguments
- `mz_input::Float64`: Target m/z value
- `atom_pool::Dict{String, Int}`: Maximum atom counts per element
- `adduct_info::AdductInfo`: Adduct information (mass, charge, name)

# Returns
- `Vector{Compound}`: All generated formulas (before tolerance filtering)
"""
function generate_formulas(mz_input::Float64, atom_pool::Dict{String, Int}, adduct_info::AdductInfo)
    formulas = Compound[]  # Type-stable array initialization
    adduct_mass = adduct_info.mass
    charge = adduct_info.charge
    
    function generate_combinations(elements, counts, idx, current)
        if idx > length(elements)
            # Extract the counts of C, H, N, and O from 'current'
            c_idx = findfirst(isequal("C"), elements)
            h_idx = findfirst(isequal("H"), elements)
            n_idx = findfirst(isequal("N"), elements)
            o_idx = findfirst(isequal("O"), elements)

            C = isnothing(c_idx) ? 0 : current[c_idx]
            H = isnothing(h_idx) ? 0 : current[h_idx]
            N = isnothing(n_idx) ? 0 : current[n_idx]
            O = isnothing(o_idx) ? 0 : current[o_idx]

            # Apply chemical heuristic rules to prune invalid formulas
            # Rule 1: Hydrogen ratio - H should be at least 0.5*C but not more than 2*C+2+N
            if C > 0 && (H < 0.5 * C || H > 2 * C + 2 + N)
                return
            end

            # Rule 2: O/C ratio - typically <= 3 for organic compounds
            if C > 0 && O > 3 * C
                return
            end

            # Rule 3: N/C ratio - typically <= 4 for organic compounds
            if C > 0 && N > 4 * C
                return
            end

            # Rule 4: DBE (Double Bond Equivalents) must be >= 0
            # DBE = (2C + 2 + N - H) / 2
            # For this to be valid: 2C + 2 + N >= H
            if 2 * C + 2 + N < H
                return
            end

            # Calculate mass and ppm
            M = sum(current[i] * ELEMENTS[elements[i]]["Relative Atomic Mass"][1] for i in eachindex(elements)) + adduct_mass - charge*0.0005485
            mz_calculated = charge == 0 ? M : M / abs(charge)

            # Construct formula string in Hill notation order
            formula = build_hill_notation(elements, current)
            ppm = ((mz_input - mz_calculated) / mz_calculated) * 1e6
            push!(formulas, Compound(formula, adduct_info.name, charge, mz_calculated, ppm))

            return
        end
        for count in 0:counts[idx]
            current[idx] = count
            generate_combinations(elements, counts, idx + 1, current)
        end
    end
    
    elements = collect(keys(atom_pool))
    n_elements = length(elements)
    max_counts = Vector{Int}(undef, n_elements)
    for i in 1:n_elements
        elem = elements[i]
        max_counts[i] = min(atom_pool[elem], round(Int, mz_input / ELEMENTS[elem]["Relative Atomic Mass"][1]))
    end
    generate_combinations(elements, max_counts, 1, zeros(Int, n_elements))
    
    return formulas
end

"""
    filter_formulas(formulas, mz_input, tolerance) -> Vector{Compound}

Filter formulas by mass accuracy tolerance.

# Arguments
- `formulas`: Vector of Compound objects to filter
- `mz_input`: Target m/z value (not currently used, kept for API consistency)
- `tolerance`: Maximum acceptable ppm error

# Returns
- `Vector{Compound}`: Formulas with |ppm| ≤ tolerance
"""
function filter_formulas(formulas, mz_input, tolerance)
    return filter(c -> abs(c.ppm) <= tolerance, formulas)
end

"""
    find_formula(mz_input::Float64; tolerance_ppm=100, atom_pool=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct="")

Find possible molecular formulas matching an experimental m/z value.

This function performs a combinatorial search through possible elemental compositions
to find formulas that match the input mass within a specified tolerance. Results are
sorted by mass accuracy (ppm error).

# Arguments
- `mz_input::Float64`: Experimental m/z value to match
- `tolerance_ppm::Number=100`: Mass tolerance in parts per million (default: 100 ppm)
- `atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10)`:
  Dictionary specifying maximum atom counts for each element to consider
- `adduct::String=""`: Adduct type - one of "H+", "Na+", "K+", "H-", "+", or "" for neutral
  - "H+": Protonated [M+H]+ (charge +1)
  - "Na+": Sodium adduct [M+Na]+ (charge +1)
  - "K+": Potassium adduct [M+K]+ (charge +1)
  - "H-": Deprotonated [M-H]- (charge -1)
  - "+": Radical cation [M]+• (charge +1, electron ionization)
  - "": Neutral molecule (charge 0)

# Returns
- `Vector{Compound}`: Array of matching formulas sorted by ppm error (best matches first)

# Examples
```julia
julia> find_formula(58.0419)
Matching formulas:
 C3H6O  m/z: 58.041865  ppm: 0.61

julia> find_formula(59.0491; adduct="H+")
Matching formulas:
 C3H6O  [M+H]+  m/z: 59.049141  ppm: -0.7

julia> find_formula(81.0284; adduct="Na+")
Matching formulas:
 C3H6O  [M+Na]+  m/z: 81.028084  ppm: 3.9

julia> find_formula(58.0419; adduct="+")
Matching formulas:
 C3H6O  [M]+  m/z: 58.041316  ppm: 1.01

julia> # Custom atom pool for peptide search
julia> find_formula(150.05; adduct="H+", atom_pool=Dict("C"=>10, "H"=>20, "N"=>5, "O"=>5, "S"=>2))
```

# Notes
- Uses chemical heuristics to prune search space (e.g., H ≥ 0.5×C)
- Larger atom pools increase computation time exponentially
- For high-resolution mass spec data, use smaller tolerance (e.g., 5-10 ppm)
- Results are printed to console and also returned as a vector

# Algorithm
1. Generates all possible formula combinations within atom_pool limits
2. Calculates theoretical m/z for each formula (including adduct/charge)
3. Filters by ppm tolerance
4. Sorts by absolute ppm error (lowest error first)

See also: [`isotopic_pattern`](@ref), [`monoisotopic_mass`](@ref), [`Compound`](@ref)
"""
function find_formula(mz_input::Float64; tolerance_ppm::Number=100, atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct::String="")
    adduct_info = get_adduct_info(adduct)
    formulas = generate_formulas(mz_input, atom_pool, adduct_info)
    matching_formulas = filter_formulas(formulas, mz_input, tolerance_ppm)

    # Sort the matching formulas by ppm
    sort!(matching_formulas, by = c -> abs(c.ppm))

    println("Matching formulas:")
    for compound in matching_formulas
        charge_sign = adduct_info.charge > 0 ? "+" : (adduct_info.charge < 0 ? "-" : "")
        adduct_string = adduct_info.name == "" ? charge_sign : "\t[" * adduct_info.name * "]" * charge_sign
        println(" ", compound.formula, adduct_string, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
    end
    return matching_formulas
end

export Compound, find_formula

end # module FindFormulaModule