"""
    Compound

Data structure representing a molecular formula match from mass spectrometry data.

# Fields
- `formula::String`: Molecular formula (e.g., "C3H6O")
- `adduct::String`: Adduct type (e.g., "M+H", "M+Na", or "" for none)
- `charge::Int`: Charge state (positive, negative, or 0 for neutral)
- `mz::Float64`: Calculated m/z value in atomic mass units
- `ppm::Float64`: Mass accuracy in parts per million (ppm)

# Examples
```julia
julia> compound = Compound("C3H6O", "M+H", 1, 59.049141, -0.7)
Compound("C3H6O", "M+H", 1, 59.049141, -0.7)
```

See also: [`findFormula`](@ref)
"""
struct Compound
    formula::String
    adduct::String
    charge::Int
    mz::Float64
    ppm::Float64
end

"""
Dictionary of common mass spectrometry adduct masses in atomic mass units (amu).

# Supported Adducts
- "M+H": Proton adduct (+1.00783 amu)
- "M+Na": Sodium adduct (+22.98977 amu)
- "M+K": Potassium adduct (+38.96371 amu)
- "M-H": Deprotonated (-1.00783 amu)
"""
const ADDUCTS = Dict("M+H" => 1.00782503223, "M+Na" => 22.989769282, "M+K" => 38.9637064864, "M-H" => -1.00782503223)

"""
    generate_formulas(mz_input, atom_pool, adduct, charge) -> Vector{Compound}

Generate all possible molecular formulas within the atom pool constraints.

Uses combinatorial search with chemical heuristics to enumerate formulas.
Currently applies H ≥ 0.5×C rule to reduce search space.

# Arguments
- `mz_input::Float64`: Target m/z value
- `atom_pool::Dict{String, Int}`: Maximum atom counts per element
- `adduct::String`: Adduct type (from ADDUCTS dictionary)
- `charge::Int`: Charge state

# Returns
- `Vector{Compound}`: All generated formulas (before tolerance filtering)
"""
function generate_formulas(mz_input::Float64, atom_pool::Dict{String, Int}, adduct::String, charge::Int)
    formulas = Compound[]  # Type-stable array initialization
    adduct_mass = get(ADDUCTS, adduct, 0.0)
    
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
            formula = join([string(elements[i], current[i] == 1 ? "" : current[i]) for i in 1:length(elements) if current[i] > 0])
            ppm = ((mz_input - mz_calculated) / mz_calculated) * 1e6
            push!(formulas, Compound(formula, adduct, charge, mz_calculated, ppm))

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
    findFormula(mz_input::Float64; tolerance_ppm=100, atom_pool=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct="", charge=0)

Find possible molecular formulas matching an experimental m/z value.

This function performs a combinatorial search through possible elemental compositions
to find formulas that match the input mass within a specified tolerance. Results are
sorted by mass accuracy (ppm error).

# Arguments
- `mz_input::Float64`: Experimental m/z value to match
- `tolerance_ppm::Number=100`: Mass tolerance in parts per million (default: 100 ppm)
- `atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10)`:
  Dictionary specifying maximum atom counts for each element to consider
- `adduct::String=""`: Adduct type - one of "M+H", "M+Na", "M+K", "M-H", or "" for none
- `charge::Int=0`: Charge state (positive or negative integer, or 0 for neutral)

# Returns
- `Vector{Compound}`: Array of matching formulas sorted by ppm error (best matches first)

# Examples
```julia
julia> findFormula(58.0419)
Matching formulas:
 C3H6O  m/z: 58.041865  ppm: 0.61

julia> findFormula(46.000; tolerance_ppm=2000, charge=1)
Matching formulas:
 CH2O2+ m/z: 46.004931  ppm: -107.18
 NO2+   m/z: 45.992355  ppm: 166.23
 ...

julia> findFormula(59.0491; adduct="M+H", charge=1)
Matching formulas:
 C3H6O  [M+H]+  m/z: 59.049141  ppm: -0.7

julia> # Custom atom pool for peptide search
julia> findFormula(150.05, atom_pool=Dict("C"=>10, "H"=>20, "N"=>5, "O"=>5, "S"=>2))
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

See also: [`isotopicPattern`](@ref), [`monoisotopicMass`](@ref), [`Compound`](@ref)
"""
function findFormula(mz_input::Float64; tolerance_ppm::Number=100, atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10), adduct::String="", charge::Int=0)
    formulas = generate_formulas(mz_input, atom_pool, adduct, charge)
    matching_formulas = filter_formulas(formulas, mz_input, tolerance_ppm)
    
    # Sort the matching formulas by ppm
    sort!(matching_formulas, by = c -> abs(c.ppm))
    
    println("Matching formulas:")
    for compound in matching_formulas
        charge_sign = charge > 0 ? "+" : (charge < 0 ? "-" : "")
        adduct_string = adduct == "" ? charge_sign : "\t[" * adduct * "]" * charge_sign
        println(" ", compound.formula, adduct_string, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
    end
    return matching_formulas
end