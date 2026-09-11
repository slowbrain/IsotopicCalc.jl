module IsotopicPatternModule

using Printf

# Access ELEMENTS from parent module
const ELEMENTS = if isdefined(parentmodule(@__MODULE__), :ELEMENTS)
    parentmodule(@__MODULE__).ELEMENTS
else
    error("ELEMENTS constant not found in parent module. Ensure IsotopicCalc is properly initialized with elements.json data.")
end

# Pre-compiled regex patterns for better performance
const ISOTOPE_TOKEN_REGEX = r"^[0-9]+[A-Z][a-z]?$"

# Formulas are ASCII tokens even when supplied through another AbstractString.
# Keep groups as count dictionaries: (CH2)1000000 needs no expanded text.
is_formula_digit(c::Char) = '0' <= c <= '9'
is_formula_upper(c::Char) = 'A' <= c <= 'Z'
is_formula_lower(c::Char) = 'a' <= c <= 'z'

function read_atom_count(formula::String, i::Int)
    start = i
    while i <= lastindex(formula) && is_formula_digit(formula[i])
        i = nextind(formula, i)
    end
    start == i && return 1, i
    count = tryparse(Int, SubString(formula, start, prevind(formula, i)))
    (isnothing(count) || count == 0) && throw(ArgumentError("Atom and group counts must be positive integers that fit in Int."))
    return count, i
end

function add_atom_count!(counts, atom, count, multiplier=1)
    try
        amount = Base.Checked.checked_mul(count, multiplier)
        counts[atom] = Base.Checked.checked_add(get(counts, atom, 0), amount)
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("Atom count for '$atom' exceeds the supported Int range."))
    end
end

"""
    parse_formula(formula::AbstractString) -> Dict{String, Int}

Read element, isotope and parenthesis tokens without expanding repeated groups.
Isotope keys omit brackets internally (e.g. "13C"); `sum_formula` restores them.
Reject malformed token boundaries, empty/unbalanced groups and count overflow.
"""
function parse_formula(input::AbstractString)
    formula = String(input)
    isempty(formula) && throw(ArgumentError("Formula cannot be empty."))
    groups = [Dict{String, Int}()]
    i = firstindex(formula)
    while i <= lastindex(formula)
        c = formula[i]
        if c == '('
            push!(groups, Dict{String, Int}())
            i = nextind(formula, i)
            continue
        elseif c == ')'
            length(groups) == 1 && throw(ArgumentError("Unmatched ')' in formula '$formula'."))
            group = pop!(groups)
            isempty(group) && throw(ArgumentError("Empty parentheses group in formula '$formula'."))
            multiplier, i = read_atom_count(formula, nextind(formula, i))
            for (atom, count) in group
                add_atom_count!(last(groups), atom, count, multiplier)
            end
            continue
        elseif c == '['
            closing = findnext(']', formula, nextind(formula, i))
            isnothing(closing) && throw(ArgumentError("Unmatched '[' in formula '$formula'."))
            atom = String(SubString(formula, nextind(formula, i), prevind(formula, closing)))
            occursin(ISOTOPE_TOKEN_REGEX, atom) || throw(ArgumentError("Malformed isotope token '[$atom]'. Use notation like [13C] or [2H]."))
            i = nextind(formula, closing)
        elseif is_formula_upper(c)
            start = i
            i = nextind(formula, i)
            if i <= lastindex(formula) && is_formula_lower(formula[i])
                i = nextind(formula, i)
            end
            atom = String(SubString(formula, start, prevind(formula, i)))
        else
            throw(ArgumentError("Unexpected character '$c' at position $i in formula '$formula'."))
        end
        count, i = read_atom_count(formula, i)
        add_atom_count!(last(groups), atom, count)
    end
    length(groups) == 1 || throw(ArgumentError("Unmatched '(' in formula '$formula'."))
    isempty(first(groups)) && throw(ArgumentError("Formula cannot be empty."))
    return first(groups)
end

# Remove only a tail whose total absolute probability fits the error budget.
# Keeping at least one state also makes normalization safe for every valid cutoff.
function prune_tail!(distribution::Dict, budget::Real)
    budget <= 0 && return distribution
    discarded = 0.0
    for (state, probability) in sort!(collect(distribution); by=last)
        length(distribution) == 1 && break
        discarded + probability > budget && break
        discarded += probability
        delete!(distribution, state)
    end
    return distribution
end

# Count-vector keys combine all permutations of an isotope composition before
# pruning. Floating-point mass sums must not identify these intermediate states.
function isotope_states(probabilities,count,budget,::Val{N}) where N
    N == 1 && return Dict(ntuple(_->count,N)=>1.0)
    states = Dict{NTuple{N,Int},Float64}(ntuple(_->0,N)=>1.0)
    for _ in 1:count
        next = Dict{NTuple{N,Int},Float64}()
        for (composition,probability) in states, i in 1:N
            probabilities[i] == 0 && continue
            key = ntuple(j->composition[j]+(i==j),N)
            next[key] = get(next,key,0.0)+probability*probabilities[i]
        end
        states = prune_tail!(next,budget)
    end
    return states
end

function element_distribution(masses,probabilities,count,budget,::Val{N}) where N
    states = isotope_states(probabilities,count,budget,Val(N))
    return [(sum(composition[i]*masses[i] for i in 1:N),probability) for (composition,probability) in states]
end

function convolve(distro1, distro2, error_budget::Real)
    distribution = Dict{Float64, Float64}()
    for (m1, a1) in distro1, (m2, a2) in distro2
        mass = m1 + m2
        distribution[mass] = get(distribution, mass, 0.0) + a1 * a2
    end
    prune_tail!(distribution, error_budget)
    return sort!([(mass, probability) for (mass, probability) in distribution])
end

function normalize_abundances(distribution::Vector{Tuple{Float64, Float64}})
    isempty(distribution) && return distribution
    max_abundance = maximum(last, distribution)
    return [(mass, abundance / max_abundance) for (mass, abundance) in distribution]
end

"""
    group_by_resolution(distribution, R)

Merge sorted peaks within `first_mass / R` of the lightest peak in each group.
Return the abundance-weighted centroid and summed abundance of each group.
Groups are anchored to their lightest peak to avoid merging a chain of peaks
whose total span exceeds the resolution interval. `R=Inf` merges equal masses only.
"""
function group_by_resolution(distribution::Vector{Tuple{Float64, Float64}}, R::Real)
    isempty(distribution) && return distribution
    sorted = sort(distribution; by=first)
    grouped = Tuple{Float64, Float64}[]
    anchor, abundance = first(sorted)
    weighted_mass = anchor * abundance
    for (mass, probability) in Iterators.drop(sorted, 1)
        if mass - anchor <= anchor / R
            weighted_mass += mass * probability
            abundance += probability
        else
            abundance > 0 && push!(grouped, (weighted_mass / abundance, abundance))
            anchor = mass
            abundance = probability
            weighted_mass = mass * probability
        end
    end
    abundance > 0 && push!(grouped, (weighted_mass / abundance, abundance))
    return grouped
end

# A single isotope lookup is shared by direct mass calculation and patterns.
function isotope_data(atom::String)
    if haskey(ELEMENTS, atom)
        data = ELEMENTS[atom]
        masses = Float64.(data["Relative Atomic Mass"])
        probabilities = Float64.(data["Isotopic Composition"])
        probabilities ./= sum(probabilities)
        return masses, probabilities
    end
    base = get_base_element(atom)
    haskey(ELEMENTS, base) || throw(ArgumentError("Unknown element '$base' (from '$atom'). Check spelling and capitalization."))
    index = findfirst(==(atom), ELEMENTS[base]["Symbol"])
    isnothing(index) && throw(ArgumentError("Unknown isotope '$atom'. Available isotopes: $(join(ELEMENTS[base]["Symbol"], ", "))."))
    return [Float64(ELEMENTS[base]["Relative Atomic Mass"][index])], [1.0]
end

function principal_mass(atom::String)::Float64
    masses, probabilities = isotope_data(atom)
    return masses[argmax(probabilities)]
end

function get_base_element(isotope::String)
  if isuppercase(isotope[end])
      return string(isotope[end])
  else
      return isotope[end-1:end]
  end
end

function sum_formula(parsed_formula::Dict{String, Int})
  elements_order = ["C", "H", "O", "N", "P", "S", "F", "Cl", "Br", "I"]  # Common elements first
  other_elements = setdiff(collect(keys(parsed_formula)), elements_order)  # Convert keys to an array
  sort!(other_elements)

  sum_formula = ""
  for element in [elements_order; other_elements]
      count = get(parsed_formula, element, 0)
      if count > 0
          sum_formula *= is_formula_digit(first(element)) ? "[$element]" : element
          sum_formula *= count > 1 ? string(count) : ""
      end
  end
  return sum_formula
end


function combine_formulas(base_formula::Dict{String, Int}, atom_delta::Dict{String, Int})
  combined_formula = deepcopy(base_formula)
  for (atom, count) in atom_delta
      new_count = Base.Checked.checked_add(get(combined_formula, atom, 0), count)
      if new_count < 0
          throw(ArgumentError("Ion transformation removes more '$atom' atoms than are present in the formula."))
      elseif new_count == 0
          delete!(combined_formula, atom)
      else
          combined_formula[atom] = new_count
      end
  end
  return combined_formula
end

include("Ions.jl")

validated_formula(formula::AbstractString) = parse_formula(formula)

include("PatternResults.jl")

"""
    isotopic_pattern(formula; abundance_cutoff=1e-5, R=10000, ion=IonSpec(),
                     annotated=false, abundances=Dict(), mixture=nothing, print=true, io=stdout)

Calculate the isotopic pattern distribution for a given chemical formula.

This function computes all possible isotopic combinations and their relative abundances
using NIST isotopic composition data. The result includes both the mass and relative
abundance of each isotopologue.

# Arguments
- `formula::AbstractString`: Chemical formula in various formats:
  - Simple: "C3H6O"
  - With parentheses: "(CH3)2CO" or "CH3COCH3"
  - With specific isotopes: "[13C]H3COCH3" or "C3[2H]6O"
  - Deuterium shorthand: "C3D6O" (equivalent to C3[2H]6O)
- `abundance_cutoff::Real=1e-5`: Final relative abundance threshold in [0, 1]; zero disables pruning
- `R::Real=10000`: Resolving power for grouping nearby peaks; `Inf` keeps distinct masses separate
- `ion=IonSpec()`: An `IonSpec` or full notation such as `[M+H]+`, `[M-H]-`, or `[M-H]+`; default is neutral
- `annotated::Bool=false`: Return PatternResult with contributors, absolute probabilities and settings
- `abundances=Dict()`: Per-element isotope probability overrides; omitted isotopes have zero probability
- `mixture=nothing`: Weighted IsotopePopulation vector with matching elemental counts
- `print::Bool=true`: Whether to print formatted output
- `io::IO=stdout`: Output destination; ignored when `print=false`

# Returns
- `Vector{Tuple{Float64, Float64}}`: Array of (mass, abundance) tuples where:
  - mass: Isotopic mass, or m/z when charged
  - abundance: Relative abundance (normalized to most abundant peak = 1.0)
- With `annotated=true`, a `PatternResult` containing peaks, probability accounting,
  calculation settings and the bundled isotope-table fingerprint is returned instead.

# Examples
```julia
julia> isotopic_pattern("CH3COCH3")
Formula:  C3H6O
----------------------------
Mass [amu]     Abundance [%]
----------------------------
58.04186481    100
59.04528960    3.351819
...

julia> isotopic_pattern("C6H12O6"; abundance_cutoff=1e-3, ion="[M+Na]+")
# Returns distribution for sodium adduct of glucose with 0.1% cutoff

julia> pattern = isotopic_pattern("C3H6O"; print=false)
# Returns data without printing (useful for programmatic access)
```

# Notes
- Uses convolution algorithm to combine isotopic distributions
- Groups span at most first_mass/R; masses are unrounded abundance-weighted centroids
- The cutoff is applied after grouping and normalization
- Internal pruning omits at most min(1e-12, abundance_cutoff × 1e-8) absolute probability
  in total, apart from floating-point roundoff; zero cutoff computes without pruning
- Explicit isotope labels stay fixed under abundance overrides
- Mixture weights and each abundance override vector must sum to one
- Annotation uses additional memory to retain isotope contributors
- Electron mass adjustments applied for charged species
- Based on NIST Atomic Weights and Isotopic Compositions database

# Throws
- `ArgumentError`: If formula contains invalid characters, brackets are unbalanced,
  or if abundance_cutoff is nonfinite/outside [0, 1] or R is NaN/nonpositive
- `ArgumentError`: If unknown element or isotope is specified

See also: [`monoisotopic_mass`](@ref), [`isotopic_pattern_protonated`](@ref), [`find_formula`](@ref)
"""
function isotopic_pattern(formula::AbstractString; abundance_cutoff=1e-5, R=10000, ion=IonSpec(),
                          annotated::Bool=false, abundances=Dict(), mixture=nothing, print::Bool=true, io::IO=stdout)
  abundance_cutoff isa Real && isfinite(abundance_cutoff) && 0 <= abundance_cutoff <= 1 || throw(ArgumentError("abundance_cutoff must be in [0,1]."))
  R isa Real && !isnan(R) && R > 0 || throw(ArgumentError("R must be positive (Inf is allowed)."))
  spec = resolve_ion(ion)
  overrides = validate_abundances(abundances)
  result = calculate_pattern(formula,spec,abundance_cutoff,R,overrides,mixture,annotated)
  formula_display = sum_formula(parse_formula(formula))
  charge = spec.charge
  final_distribution = annotated ? [(p.mass,p.abundance) for p in result.peaks] : result
  if print
    # Format the output
    output = "\nFormula: \e[1;32m $formula_display"
    notation = string(spec)
    notation == "[M]" || (output *= ".$notation")
    output *= "\n\e[m----------------------------\n"
    # constant padding in title independent on tab length
    pad_tileIO = IOBuffer(); @printf(pad_tileIO, "%5s", ""); pad_title = String(take!(pad_tileIO))
    output *= (charge == 0 ? "Mass [amu]" : "m/z")*pad_title*"Abundance [%]\n"
    output *= "----------------------------\n"
    for (m, a) in final_distribution
        mm = @sprintf("%.8f", m)
        aa = @sprintf("%.8g", a * 100)

        # padding list dependent on mass
        format_pad = Printf.Format("%"*string(15-length(mm)+length(aa))*"s")
        padIO = IOBuffer(); Printf.format(padIO, format_pad, aa); pad = String(take!(padIO))

        output *= mm * pad *"\n"
    end
    output *= "Found $(length(final_distribution)) isotopic masses for $abundance_cutoff abundance limit."

    println(io, output)
  end
  return result
end

# special cases
"""
    isotopic_pattern_protonated(formula::AbstractString; kwargs...)

Convenience function to calculate the isotopic pattern for a protonated molecule [M+H]⁺.

This is equivalent to calling `isotopic_pattern(formula; ion="[M+H]+", kwargs...)`.
Keywords such as `R`, `abundance_cutoff`, `print`, and `io` are forwarded.

# Arguments
- `formula::AbstractString`: Chemical formula (same format as `isotopic_pattern`)

# Returns
- `Vector{Tuple{Float64, Float64}}`: Array of (mass, abundance) tuples for [M+H]⁺

# Examples
```julia
julia> isotopic_pattern_protonated("C3H6O")
# Calculates pattern for protonated acetone (m/z 59.049)
```

See also: [`isotopic_pattern`](@ref), [`monoisotopic_mass_protonated`](@ref)
"""
function isotopic_pattern_protonated(formula::AbstractString; kwargs...)
    haskey(kwargs, :ion) && throw(ArgumentError("The protonated convenience function fixes ion=\"[M+H]+\"."))
    return isotopic_pattern(formula; ion="[M+H]+", kwargs...)
end

"""
    monoisotopic_mass(formula::AbstractString)

Calculate the monoisotopic mass (most abundant isotopic composition) for a given formula.

Sums the mass of the most abundant isotope of each element (or the specified isotope),
without calculating the full isotopic distribution or printing output.

# Arguments
- `formula::AbstractString`: Chemical formula (same format as `isotopic_pattern`)

# Returns
- `Float64`: Monoisotopic mass in atomic mass units (amu)

# Examples
```julia
julia> monoisotopic_mass("CH3COCH3")
58.041864812949996

julia> monoisotopic_mass("C6H12O6")
180.06338810418

julia> monoisotopic_mass("[13C]H3COCH3")
59.04521964802  # One carbon is ¹³C instead of ¹²C
```

# Notes
- Computes directly from atom counts without generating an isotope pattern
- Returns unrounded exact masses, independent of resolution and abundance cutoffs
- Uses natural isotopic abundances from NIST database
- For charged species, use the `ion` parameter in `isotopic_pattern` instead

See also: [`isotopic_pattern`](@ref), [`monoisotopic_mass_protonated`](@ref)
"""
function monoisotopic_mass(formula::AbstractString)
    parsed = validated_formula(formula)
    return sum(parsed[atom] * principal_mass(atom) for atom in sort!(collect(keys(parsed))))
end

"""
    monoisotopic_mass_protonated(formula::AbstractString)

Calculate the monoisotopic mass of a protonated molecule [M+H]⁺.

Adds the principal hydrogen isotope mass and subtracts the electron mass directly.
This does not depend on spectral resolution, peak merging, or abundance cutoffs.
Useful for mass spectrometry applications where positive-mode ESI is common.

# Arguments
- `formula::AbstractString`: Chemical formula (same format as `isotopic_pattern`)

# Returns
- `Float64`: Monoisotopic m/z value for [M+H]⁺ in atomic mass units

# Examples
```julia
julia> monoisotopic_mass_protonated("C3H6O")
59.04914134518  # Acetone + proton

julia> monoisotopic_mass_protonated("C6H12O6")
181.07066463641  # Glucose + proton
```

See also: [`monoisotopic_mass`](@ref), [`isotopic_pattern_protonated`](@ref)
"""
monoisotopic_mass_protonated(formula::AbstractString) = monoisotopic_mass(formula) + principal_mass("H") - principal_mass("e")

export isotopic_pattern, monoisotopic_mass, isotopic_pattern_protonated, monoisotopic_mass_protonated

end # module IsotopicPatternModule
