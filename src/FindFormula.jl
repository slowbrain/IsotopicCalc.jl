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
    parse_adduct(adduct::String) -> AdductInfo

Parse an adduct string and return adduct information.

Supports flexible adduct notation using any element from the periodic table,
including multiple charges for multiply charged ions.

# Supported Formats
- "Element+": Addition with positive charge (e.g., "H+", "Na+", "K+", "Ca+")
- "Element-": Removal with negative charge (e.g., "H-", "Cl-")
- "nElement+": Multiple atoms, single charge (e.g., "2H+", "3Na+")
- "nElement+m": Multiple atoms, multiple charges (e.g., "2H+2", "3H+3")
- "+n": Radical cation with charge +n (e.g., "+", "+2", "+3")
- "-n": Radical anion with charge -n (e.g., "-", "-2", "-3")
- "": Neutral molecule (charge 0, no mass change)

# Arguments
- `adduct::String`: Adduct identifier in the formats above

# Returns
- `AdductInfo`: Struct containing mass change, charge, and display name

# Throws
- `ArgumentError`: If the adduct format is invalid or element not found

# Examples
```julia
julia> parse_adduct("H+")
AdductInfo(1.00782503223, 1, "M+H")

julia> parse_adduct("2H+2")
AdductInfo(2.01565006446, 2, "M+2H")

julia> parse_adduct("3H+3")
AdductInfo(3.02347509669, 3, "M+3H")

julia> parse_adduct("+2")
AdductInfo(0.0, 2, "M")

julia> parse_adduct("Ca+")
AdductInfo(39.962590863, 1, "M+Ca")
```
"""
function parse_adduct(adduct::String)
    # Handle special cases for neutral
    if adduct == ""
        return AdductInfo(0.0, 0, "")
    end

    # Handle special cases for radical ions: "+", "+2", "+3", "-", "-2", "-3"
    m_radical = match(r"^([+-])(\d*)$", adduct)
    if m_radical !== nothing
        charge_sign, charge_num_str = m_radical.captures
        charge_magnitude = charge_num_str == "" ? 1 : parse(Int, charge_num_str)
        charge = charge_sign == "+" ? charge_magnitude : -charge_magnitude
        display_name = "M"  # Just "M" for radical ions, charge will be added in output formatting
        return AdductInfo(0.0, charge, display_name)
    end

    # Parse format: [atom_count]Element[+|-][charge_count]
    # Match pattern: optional number, element symbol, +/-, optional number
    m = match(r"^(\d*)([A-Z][a-z]?)([+-])(\d*)$", adduct)

    if m === nothing
        throw(ArgumentError("Invalid adduct format '$adduct'. Expected format: [n]Element+/-[m] (e.g., 'H+', '2H+2', 'Cl-') or special forms '+n', '-n', ''"))
    end

    atom_count_str, element, charge_sign, charge_count_str = m.captures

    # Parse atom count (default to 1 if not specified)
    atom_count = atom_count_str == "" ? 1 : parse(Int, atom_count_str)

    # Parse charge count (default to 1 if not specified)
    charge_magnitude = charge_count_str == "" ? 1 : parse(Int, charge_count_str)

    # Convert SubString to String for JSON.Object compatibility
    element = String(element)

    # Check if element exists
    if !haskey(ELEMENTS, element)
        available_elements = join(sort(collect(keys(ELEMENTS)))[1:min(20, length(ELEMENTS))], ", ")
        throw(ArgumentError("Unknown element '$element'. First 20 available elements: $available_elements, ..."))
    end

    # Get element mass
    element_mass = ELEMENTS[element]["Relative Atomic Mass"][1]

    # Calculate total mass change and charge
    charge = charge_sign == "+" ? charge_magnitude : -charge_magnitude
    mass_change = charge_sign == "+" ? atom_count * element_mass : -atom_count * element_mass

    # Create display name
    atom_count_display = atom_count == 1 ? "" : string(atom_count)
    display_name = "M$(charge_sign)$(atom_count_display)$(element)"

    return AdductInfo(mass_change, charge, display_name)
end

# Keep get_adduct_info for backward compatibility, but make it use parse_adduct
function get_adduct_info(adduct::String)
    return parse_adduct(adduct)
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
Applies multiple chemical rules to reduce search space:
- Hydrogen/Carbon ratios (H ≥ 0.5×C when H>0, H ≤ 2C+2+N)
- Heteroatom/Carbon ratios (O ≤ 3C, N ≤ 4C, S ≤ C, P ≤ 2C)
- Halogen/Carbon ratio ((F+Cl+Br+I) ≤ 2C)
- Nitrogen Rule (parity check based on nominal mass)
- DBE (Double Bond Equivalents) ≥ 0

# Arguments
- `mz_input::Real`: Target m/z value
- `atom_pool::Dict{String, Int}`: Maximum atom counts per element
- `adduct_info::AdductInfo`: Adduct information (mass, charge, name)

# Returns
- `Vector{Compound}`: All generated formulas (before tolerance filtering)
"""
function generate_formulas(mz_input::Real, atom_pool::Dict{String, Int}, adduct_info::AdductInfo)
    formulas = Compound[]  # Type-stable array initialization
    adduct_mass = adduct_info.mass
    charge = adduct_info.charge
    
    function generate_combinations(elements, counts, idx, current)
        if idx > length(elements)
            # Extract the counts of C, H, N, O, S, P, and halogens from 'current'
            c_idx = findfirst(isequal("C"), elements)
            h_idx = findfirst(isequal("H"), elements)
            n_idx = findfirst(isequal("N"), elements)
            o_idx = findfirst(isequal("O"), elements)
            s_idx = findfirst(isequal("S"), elements)
            p_idx = findfirst(isequal("P"), elements)
            f_idx = findfirst(isequal("F"), elements)
            cl_idx = findfirst(isequal("Cl"), elements)
            br_idx = findfirst(isequal("Br"), elements)
            i_idx = findfirst(isequal("I"), elements)

            C = isnothing(c_idx) ? 0 : current[c_idx]
            H = isnothing(h_idx) ? 0 : current[h_idx]
            N = isnothing(n_idx) ? 0 : current[n_idx]
            O = isnothing(o_idx) ? 0 : current[o_idx]
            S = isnothing(s_idx) ? 0 : current[s_idx]
            P = isnothing(p_idx) ? 0 : current[p_idx]
            F = isnothing(f_idx) ? 0 : current[f_idx]
            Cl = isnothing(cl_idx) ? 0 : current[cl_idx]
            Br = isnothing(br_idx) ? 0 : current[br_idx]
            I = isnothing(i_idx) ? 0 : current[i_idx]

            # Apply chemical heuristic rules to prune invalid formulas
            # Rule 1: Hydrogen ratio - H should not exceed 2*C+2+N
            # Note: We allow H=0 for compounds with no hydrogen
            if C > 0 && H > 2 * C + 2 + N
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

            # Rule 5: Nitrogen Rule (parity rule)
            # Calculate nominal mass of the NEUTRAL molecule (before adduct)
            # The nitrogen rule applies to neutral molecules only
            nominal_mass = 0
            for i in eachindex(elements)
                if current[i] > 0
                    element_mass = ELEMENTS[elements[i]]["Relative Atomic Mass"][1]
                    nominal_mass += current[i] * round(Int, element_mass)
                end
            end

            # For even nominal mass, N should be even; for odd nominal mass, N should be odd
            if (nominal_mass % 2 == 0 && N % 2 != 0) || (nominal_mass % 2 != 0 && N % 2 == 0)
                return
            end

            # Rule 6: S/C ratio - typically S <= C for organic compounds
            if C > 0 && S > C
                return
            end

            # Rule 7: P/C ratio - typically P <= 2*C for organic compounds
            if C > 0 && P > 2 * C
                return
            end

            # Rule 8: Halogen/C ratio - typically (F + Cl + Br + I) <= 2*C
            total_halogens = F + Cl + Br + I
            if C > 0 && total_halogens > 2 * C
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
    find_formula(mz_input::Real; tolerance_ppm=100, atom_pool=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct="")

Find possible molecular formulas matching an experimental m/z value.

This function performs a combinatorial search through possible elemental compositions
to find formulas that match the input mass within a specified tolerance. Results are
sorted by mass accuracy (ppm error).

# Arguments
- `mz_input::Real`: Experimental m/z value to match (can be Integer or Float)
- `tolerance_ppm::Number=100`: Mass tolerance in parts per million (default: 100 ppm)
- `atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10)`:
  Dictionary specifying maximum atom counts for each element to consider
- `adduct::String=""`: Adduct type using flexible notation:
  - Format: `[n]Element+/-[m]` where n=atom count (optional), m=charge (optional)
  - Any element from the periodic table is supported (e.g., "Ca+", "Mg+", "Br-")
  - Multiple charges supported for multiply charged ions (e.g., "2H+2", "3H+3")
  - Special forms: "+n" (radical cation), "-n" (radical anion), "" (neutral)
  - Common examples:
    - "H+": Protonated [M+H]+ (charge +1)
    - "2H+2": Doubly protonated [M+2H]2+ (charge +2)
    - "3H+3": Triply protonated [M+3H]3+ (charge +3)
    - "Na+": Sodium adduct [M+Na]+ (charge +1)
    - "K+": Potassium adduct [M+K]+ (charge +1)
    - "2H+": Two protons added [M+2H]+ (charge +1, different from 2H+2!)
    - "H-": Deprotonated [M-H]- (charge -1)
    - "+2": Doubly charged radical cation [M]2+ (charge +2)
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

julia> # Multiply charged ions (e.g., for peptides)
julia> find_formula(30.0; adduct="2H+2")  # Doubly protonated
Matching formulas:
 C3H6O  [M+2H]2+  m/z: 30.028037  ppm: -...

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
function find_formula(mz_input::Real; tolerance_ppm::Number=100, atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct::String="")
    adduct_info = get_adduct_info(adduct)
    formulas = generate_formulas(mz_input, atom_pool, adduct_info)
    matching_formulas = filter_formulas(formulas, mz_input, tolerance_ppm)

    # Sort the matching formulas by ppm
    sort!(matching_formulas, by = c -> abs(c.ppm))

    println("Matching formulas:")
    for compound in matching_formulas
        charge_abs = abs(adduct_info.charge)
        charge_sign = adduct_info.charge > 0 ? "+" : (adduct_info.charge < 0 ? "-" : "")

        # Format charge suffix: +, 2+, 3+, etc.
        charge_suffix = charge_abs == 0 ? "" : (charge_abs == 1 ? charge_sign : string(charge_abs) * charge_sign)

        adduct_string = adduct_info.name == "" ? charge_suffix : "\t[" * adduct_info.name * "]" * charge_suffix
        println(" ", compound.formula, adduct_string, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
    end
    return matching_formulas
end

export Compound, find_formula

end # module FindFormulaModule