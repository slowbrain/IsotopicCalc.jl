module FindFormulaModule

using ..IsotopicPatternModule: principal_mass, IonSpec, resolve_ion, ion_delta

# Access ELEMENTS from parent module
const ELEMENTS = if isdefined(parentmodule(@__MODULE__), :ELEMENTS)
    parentmodule(@__MODULE__).ELEMENTS
else
    error("ELEMENTS constant not found in parent module. Ensure IsotopicCalc is properly initialized with elements.json data.")
end

"""
    Compound

A molecular formula match with `formula` in Hill notation, full `ion` notation
(e.g. `[M+H]+` or `[M-H]-`), signed `charge`, calculated `mz`, and mass error
`ppm = 1e6 * (observed_mz - calculated_mz) / calculated_mz`.

See also: [`find_formula`](@ref), [`IonSpec`](@ref)
"""
struct Compound
    formula::String
    ion::String
    charge::Int
    mz::Float64
    ppm::Float64
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

include("SearchResults.jl")

"""
    generate_formulas(mz_input, atom_pool, spec; heuristics=true, tolerance_ppm=0) -> Vector{Compound}

Generate matching molecular formulas using neutral-mass bounds during enumeration.

Uses combinatorial search with chemical heuristics to enumerate formulas.
Applies integral, nonnegative DBE to closed-shell CHNO/halogen compositions,
counting halogens and deuterium as hydrogen equivalents. Other elements bypass
this restricted valence model. With `heuristics=true`, also applies empirical
ratios O ≤ 3C, N ≤ 4C, S ≤ C, P ≤ 2C when carbon is present.
Atom removal is always checked against candidate atom counts.

# Arguments
- `mz_input::Real`: Target m/z value
- `atom_pool::Dict{String, Int}`: Maximum atom counts per element
- `spec::IonSpec`: Atom additions, losses, and signed charge

# Returns
- `Vector{Compound}`: Matching formulas within tolerance (before sorting)
"""
function generate_formulas(mz_input::Real, atom_pool::Dict{String, Int}, spec::IonSpec; heuristics::Bool=true, tolerance_ppm::Real=0, tolerance_da::Real=0,
                           constraints::SearchConstraints=SearchConstraints(), mode::Symbol=:neutral, stats::SearchDiagnostics=SearchDiagnostics())
    formulas = Compound[]
    delta = ion_delta(spec)
    atom_mass_change = sum(Float64[delta[atom]*principal_mass(atom) for atom in sort!(collect(keys(delta)))])
    charge = spec.charge
    notation = string(spec)
    electron_mass = principal_mass("e")
    charge_scale = max(abs(charge), 1)

    # |(observed - theoretical) / theoretical| <= tolerance, for positive m/z.
    fraction = tolerance_ppm / 1e6
    lower_mass = charge_scale * min(mz_input / (1 + fraction),mz_input-tolerance_da) - atom_mass_change + charge * electron_mass
    upper_mass = fraction < 1 ? charge_scale * max(mz_input / (1 - fraction),mz_input+tolerance_da) - atom_mass_change + charge * electron_mass : Inf
    # Branch arithmetic can add masses in a different order than the final ppm
    # calculation. Widen only the pruning bounds; final acceptance stays exact.
    scale = max(abs(lower_mass), isfinite(upper_mass) ? abs(upper_mass) : 0.0, abs(atom_mass_change), 1.0)
    slack = 64 * eps(Float64) * scale
    lower_mass -= slack
    upper_mass += slack

    data = Tuple{String, Int, Float64}[]
    for (element, limit) in atom_pool
        limit >= 0 || throw(ArgumentError("atom_pool counts must be non-negative (got $limit for '$element')."))
        haskey(ELEMENTS, element) && element != "e" || throw(ArgumentError("Unknown atom_pool element '$element'."))
        push!(data, (element, limit, principal_mass(element)))
    end
    # Heavy elements first leave a narrow feasible count interval for hydrogen.
    sort!(data; by=item -> (-item[3], item[1]))
    elements = [item[1] for item in data]
    masses = [item[3] for item in data]
    n = length(elements)
    (n == 0 || upper_mass <= 0) && return formulas
    minimum_counts = [max(0, -get(delta, element, 0)) for element in elements]
    maximum_counts = [item[2] for item in data]
    for (element, count_change) in delta
        if count_change < 0 && !haskey(atom_pool,element)
            stats.status = :infeasible_ion
            return formulas
        end
    end
    for i in 1:n
        bound = max(0.0, upper_mass / masses[i])
        if bound < maximum_counts[i]
            maximum_counts[i] = floor(Int, bound)
        end
        if minimum_counts[i] > maximum_counts[i]
            stats.status = :infeasible_ion
            return formulas
        end
    end

    suffix_min = zeros(Float64, n+1)
    suffix_max = zeros(Float64, n+1)
    for i in n:-1:1
        suffix_min[i] = suffix_min[i+1] + minimum_counts[i] * masses[i]
        suffix_max[i] = suffix_max[i+1] + maximum_counts[i] * masses[i]
    end
    current = zeros(Int, n)
    # Compute final masses in canonical symbol order, independent of search order.
    mass_order = sortperm(elements)
    indices = Dict(element => i for (i, element) in enumerate(elements))
    c_idx, h_idx, d_idx, n_idx, o_idx, s_idx, p_idx, f_idx, cl_idx, br_idx, i_idx =
        (get(indices, element, 0) for element in ("C", "H", "D", "N", "O", "S", "P", "F", "Cl", "Br", "I"))
    restricted_valence = [element in ("C", "H", "D", "N", "O", "F", "Cl", "Br", "I") for element in elements]
    atom_count(index) = index == 0 ? 0 : current[index]

    function visit(idx, partial_mass)
        if idx > n
            all(iszero, current) && return
            neutral_mass = sum(current[i] * masses[i] for i in mass_order)
            mz = (neutral_mass + atom_mass_change - charge * electron_mass) / charge_scale
            mz > 0 || return
            ppm = ((mz_input - mz) / mz) * 1e6
            (abs(ppm) <= tolerance_ppm || abs(mz_input-mz) <= tolerance_da) || return
            stats.mass_matches += 1
            for (element,minimum) in constraints.min_atoms
                atom_count(get(indices,element,0)) < minimum && return reject!(stats,:minimum_atoms)
            end

            C, H, D, N, O, S, P = (atom_count(i) for i in (c_idx, h_idx, d_idx, n_idx, o_idx, s_idx, p_idx))
            supported = all(i->current[i]==0 || restricted_valence[i],1:n)
            halogens = atom_count(f_idx)+atom_count(cl_idx)+atom_count(br_idx)+atom_count(i_idx)
            twice_dbe = 2*C+2+N-H-D-halogens+(mode == :ion ? charge : 0)
            if supported && constraints.valence == :auto
                (twice_dbe < 0 || isodd(twice_dbe)) && return reject!(stats,:valence)
            end
            if !isnothing(constraints.dbe)
                supported || return reject!(stats,:unsupported_dbe)
                constraints.dbe[1] <= twice_dbe/2 <= constraints.dbe[2] || return reject!(stats,:dbe)
            end
            for ((numerator,denominator),(lo,hi)) in constraints.ratios
                den = atom_count(get(indices,denominator,0))
                den == 0 && continue
                ratio = atom_count(get(indices,numerator,0))/den
                lo <= ratio <= hi || return reject!(stats,:ratios)
            end
            if heuristics && C > 0
                for (element,count,limit) in (("O",O,3),("N",N,4),("S",S,1),("P",P,2))
                    !haskey(constraints.ratios,(element,"C")) && count > limit*C && return reject!(stats,:heuristics)
                end
            end
            formula = build_hill_notation(elements, current)
            push!(formulas, Compound(formula, notation, charge, mz, ppm))
            return
        end
        # Remaining elements bound what each branch can still reach. Intersect
        # that interval with the allowed count range before enumerating it.
        low = (lower_mass - partial_mass - suffix_max[idx+1]) / masses[idx]
        high = (upper_mass - partial_mass - suffix_min[idx+1]) / masses[idx]
        (low > maximum_counts[idx] || high < minimum_counts[idx]) && return
        first_count = low <= minimum_counts[idx] ? minimum_counts[idx] : ceil(Int, low)
        last_count = high >= maximum_counts[idx] ? maximum_counts[idx] : floor(Int, high)
        for count in first_count:last_count
            current[idx] = count
            visit(idx+1, partial_mass + count * masses[idx])
        end
    end
    visit(1, 0.0)
    stats.status = !isempty(formulas) ? :matched : stats.mass_matches > 0 ? :filtered_out : :no_mass_match
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
    find_formula(mz_input; tolerance_ppm=100, tolerance_da=0, atom_pool=...,
                 ion=IonSpec(), mode=:neutral, constraints=SearchConstraints(),
                 diagnostics=false, heuristics=true, print=true, io=stdout)

Find possible molecular formulas matching an experimental m/z value.

This function performs a combinatorial search through possible elemental compositions
to find formulas that match the input mass within a specified tolerance. Results are
sorted by mass accuracy (ppm error).

# Arguments
- `mz_input::Real`: Experimental m/z value to match (can be Integer or Float)
- `tolerance_ppm::Real=100`: Mass tolerance in parts per million
- `tolerance_da::Real=0`: Absolute tolerance floor; acceptance uses the larger of this and ppm
- `ion=IonSpec()`: An `IonSpec` or full notation such as `[M-H]-` or `[M-H]+`; default is neutral
- `mode::Symbol=:neutral`: Neutral precursor search, or `:ion` for an actual charged composition
- `constraints::SearchConstraints`: Minimum counts, explicit ratio/DBE intervals, valence mode
- `diagnostics::Bool=false`: Return SearchResult with status, rejection counts and settings
- `atom_pool::Dict{String, Int}=Dict("C"=>20, "H"=>100, "O"=>10, "N"=>10)`:
  Dictionary specifying maximum atom counts for each element to consider
- `heuristics::Bool=true`: Apply empirical organic composition limits
- `print::Bool=true`: Whether to print matching formulas
- `io::IO=stdout`: Output destination; ignored when `print=false`

# Returns
- `Vector{Compound}`: Array of matching formulas sorted by ppm error (best matches first)
- `diagnostics=true` returns a `SearchResult` whose `.matches` contains that vector

# Examples
```julia
julia> find_formula(58.0419)
Matching formulas:
 C3H6O  m/z: 58.041865  ppm: 0.61

julia> find_formula(59.0491; ion="[M+H]+")
Matching formulas:
 C3H6O  [M+H]+  m/z: 59.049141  ppm: -0.7

julia> find_formula(81.0311; ion="[M+Na]+")
Matching formulas:
 C3H6O  [M+Na]+  m/z: 81.031086  ppm: 0.18

julia> find_formula(58.0419; ion="[M]+")
Matching formulas:
 C3H6O  [M]+  m/z: 58.041316  ppm: 10.06

julia> # Multiply charged ions (e.g., for peptides)
julia> find_formula(30.0282; ion="[M+2H]2+")  # Doubly protonated
Matching formulas:
 C3H6O  [M+H2]2+  m/z: 30.028209  ppm: -0.3

julia> # Custom atom pool for peptide search
julia> find_formula(150.05; ion="[M+H]+", atom_pool=Dict("C"=>10, "H"=>20, "N"=>5, "O"=>5, "S"=>2))
```

# Notes
- `heuristics=false` disables default empirical O/C, N/C, S/C, and P/C limits
- Explicit constraints still apply; explicit ratios replace matching default limits
- Ion mode requires a charge-only spec. Its default CHNO/halogen valence model
  assumes even-electron ions, adding charge to the DBE numerator. Use
  `SearchConstraints(valence=:none)` for radical ions or other assumptions.
- Diagnostic mass-match counts are measured after atom-pool/atom-loss availability
  bounds. They distinguish no mass matches from candidates rejected by constraints.
- Closed-shell CHNO/halogen candidates must have integral, nonnegative DBE;
  deuterium counts as hydrogen. Other element sets bypass this valence model.
- Masses use the most abundant isotope of each element without rounding
- Atom removal cannot exceed the candidate's available atom counts
- Larger atom pools increase computation time exponentially
- For high-resolution mass spec data, use smaller tolerance (e.g., 5-10 ppm)
- Results are returned as a vector; `print=false` suppresses output and `io` redirects it

# Algorithm
1. Bounds partial compositions by the neutral mass interval and remaining atom counts
2. Checks exact ppm and chemical constraints before constructing formula strings
3. Sorts by absolute ppm error (lowest error first)

See also: [`isotopic_pattern`](@ref), [`monoisotopic_mass`](@ref), [`Compound`](@ref)
"""
function find_formula(mz_input::Real; tolerance_ppm::Real=100, tolerance_da::Real=0,
                      atom_pool::Dict{String,Int}=Dict("C"=>20,"H"=>100,"O"=>10,"N"=>10),
                      ion=IonSpec(), mode::Symbol=:neutral,
                      constraints::SearchConstraints=SearchConstraints(), diagnostics::Bool=false,
                      heuristics::Bool=true, print::Bool=true, io::IO=stdout)
    isfinite(mz_input) && mz_input > 0 || throw(ArgumentError("mz_input must be finite and positive."))
    isfinite(tolerance_ppm) && tolerance_ppm >= 0 || throw(ArgumentError("tolerance_ppm must be finite and non-negative."))
    isfinite(tolerance_da) && tolerance_da >= 0 || throw(ArgumentError("tolerance_da must be finite and non-negative."))
    mode in (:neutral,:ion) || throw(ArgumentError("mode must be :neutral or :ion."))
    spec = resolve_ion(ion)
    if mode == :ion
        spec.charge != 0 && isempty(spec.additions) && isempty(spec.losses) || throw(ArgumentError("Ion-composition searches require a charge-only ion specification."))
    end
    stats = SearchDiagnostics()
    matching_formulas = generate_formulas(mz_input,atom_pool,spec;heuristics=heuristics,
        tolerance_ppm=tolerance_ppm,tolerance_da=tolerance_da,constraints=constraints,mode=mode,stats=stats)

    # Sort the matching formulas by ppm
    sort!(matching_formulas, by = c -> abs(c.ppm))

    if print
        println(io, "Matching formulas:")
        for compound in matching_formulas
            ion_display = compound.ion == "[M]" ? "" : "\t" * compound.ion
            println(io, " ", compound.formula, ion_display, "\tm/z: ", string(round(compound.mz, digits=6)), "\tppm: ", string(round(compound.ppm, digits=2)))
        end
    end
    if diagnostics
        settings = (mz_input=mz_input,tolerance_ppm=tolerance_ppm,tolerance_da=tolerance_da,
                    atom_pool=copy(atom_pool),ion=deepcopy(spec),mode=mode,constraints=deepcopy(constraints),
                    heuristics=heuristics,package_version=string(parentmodule(@__MODULE__).PACKAGE_VERSION),
                    isotope_data_sha256=parentmodule(@__MODULE__).ISOTOPE_DATA_SHA256)
        return SearchResult(matching_formulas,stats,settings)
    end
    return matching_formulas
end

export Compound, find_formula

end # module FindFormulaModule
