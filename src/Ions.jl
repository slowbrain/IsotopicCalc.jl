"""
    IonSpec(; add="", remove="", charge=0)
    IonSpec(notation)

Ion transformation: atom additions and losses, independent of signed charge.
Use full `[M+...-...]z+/-` notation, or omit the charge suffix for a neutral
transformation. `[M-H]-` removes hydrogen with charge -1; `[M-H]+` removes
hydrogen with charge +1; `[M-H]` removes hydrogen without charging.
`[M+Cl]-` adds chlorine with charge -1. `[M]` (or `IonSpec()`) changes nothing.
A charge-only spec such as `[M]2+` describes a composition already containing
every ion atom. Compact shorthand such as `H-`, `Cl-`, or `+2` is not accepted.

`string(spec)` returns full notation derived from the composition and charge.
"""
struct IonSpec
    additions::Dict{String,Int}
    losses::Dict{String,Int}
    charge::Int
    function IonSpec(additions::Dict{String,Int}, losses::Dict{String,Int}, charge::Int)
        charge == typemin(Int) && throw(ArgumentError("Unsupported charge magnitude."))
        for composition in (additions, losses), (atom,count) in composition
            count > 0 || throw(ArgumentError("Ion atom counts must be positive."))
            principal_mass(atom)
        end
        new(copy(additions),copy(losses),charge)
    end
end

function IonSpec(; add::AbstractString="", remove::AbstractString="", charge::Integer=0)
    additions = isempty(add) ? Dict{String,Int}() : parse_formula(add)
    losses = isempty(remove) ? Dict{String,Int}() : parse_formula(remove)
    return IonSpec(additions,losses,Int(charge))
end

function Base.string(ion::IonSpec)
    added = isempty(ion.additions) ? "" : "+" * sum_formula(ion.additions)
    removed = isempty(ion.losses) ? "" : "-" * sum_formula(ion.losses)
    suffix = ion.charge == 0 ? "" :
        (abs(ion.charge) == 1 ? "" : string(abs(ion.charge))) * (ion.charge > 0 ? "+" : "-")
    return "[M" * added * removed * "]" * suffix
end

function ion_term!(composition, term)
    m = match(r"^([0-9]*)(.+)\z",term)
    isnothing(m) && throw(ArgumentError("Empty ion term."))
    digits, formula = m.captures
    multiplier = isempty(digits) ? 1 : tryparse(Int,digits)
    (isnothing(multiplier) || multiplier <= 0) && throw(ArgumentError("Invalid ion multiplier."))
    for (atom,count) in parse_formula(formula)
        add_atom_count!(composition,atom,count,multiplier)
    end
end

function IonSpec(notation::AbstractString)
    text = String(notation)
    additions, losses = Dict{String,Int}(), Dict{String,Int}()
    full = match(r"^\[M(.*)\](?:([0-9]*)([+-]))?\z",text)
    isnothing(full) && throw(ArgumentError("Invalid ion notation '$text'. Use [M+H]+, [M-H]-, [M-H]+, [M+Cl]-, [M]2+, or [M] for neutral."))
    terms,digits,sign = full.captures
    position = 1
    for term in eachmatch(r"([+-])([^+-]+)",terms)
        term.offset == position || throw(ArgumentError("Invalid ion notation '$text'."))
        ion_term!(term[1] == "+" ? additions : losses,term[2])
        position += ncodeunits(term.match)
    end
    position == ncodeunits(terms)+1 || throw(ArgumentError("Invalid ion notation '$text'."))
    isnothing(sign) && return IonSpec(additions,losses,0)
    magnitude = isempty(digits) ? 1 : tryparse(Int,digits)
    (isnothing(magnitude) || magnitude <= 0) && throw(ArgumentError("Ion charge must be nonzero and fit in Int."))
    return IonSpec(additions,losses,sign == "+" ? magnitude : -magnitude)
end

function ion_delta(ion::IonSpec)
    delta = copy(ion.additions)
    for (atom,count) in ion.losses
        delta[atom] = Base.Checked.checked_sub(get(delta,atom,0),count)
        delta[atom] == 0 && delete!(delta,atom)
    end
    return delta
end

resolve_ion(ion::AbstractString) = IonSpec(ion)
resolve_ion(ion::IonSpec) = IonSpec(ion.additions,ion.losses,ion.charge)
resolve_ion(ion) = throw(ArgumentError("ion must be an IonSpec or full [M...] notation."))

function ion_composition(formula, ion::IonSpec)
    base = parse_formula(formula)
    for atom in keys(base)
        principal_mass(atom)
    end
    return combine_formulas(base,ion_delta(ion))
end

"""
    ion_mass(formula; ion=IonSpec())

Unrounded principal-isotope mass, or m/z for a charged ion. Adds/removes atoms
before applying electron mass and charge division; no spectrum is generated.
Abundance enrichment does not change these reference principal-isotope masses.
"""
function ion_mass(formula::AbstractString; ion=IonSpec())
    spec = resolve_ion(ion)
    composition = ion_composition(formula,spec)
    mass = sum(Float64[composition[a]*principal_mass(a) for a in sort!(collect(keys(composition)))])
    mass -= spec.charge*principal_mass("e")
    mass > 0 || throw(ArgumentError("Ion composition must have positive mass."))
    return mass/max(abs(spec.charge),1)
end
