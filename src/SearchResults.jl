"""
    SearchConstraints(; min_atoms=Dict(), ratios=Dict(), dbe=nothing, valence=:auto)

Minimum counts, ratio intervals `("O","C") => (0,4)`, and an optional DBE
interval. Explicit ratios replace the corresponding default heuristic. Ratios
with a zero denominator are ignored. `valence=:none` disables automatic valence
rejection; an explicit DBE interval is still applied.
"""
struct SearchConstraints
    min_atoms::Dict{String,Int}
    ratios::Dict{Tuple{String,String},Tuple{Float64,Float64}}
    dbe::Union{Nothing,Tuple{Float64,Float64}}
    valence::Symbol
    function SearchConstraints(;min_atoms=Dict(),ratios=Dict(),dbe=nothing,valence=:auto)
        valence in (:auto,:none) || throw(ArgumentError("valence must be :auto or :none."))
        minimums = Dict{String,Int}()
        for (atom,count) in min_atoms
            element = String(atom)
            haskey(ELEMENTS,element) && element != "e" || throw(ArgumentError("Unknown minimum-count element '$element'."))
            count isa Integer && count >= 0 || throw(ArgumentError("Minimum counts must be nonnegative integers."))
            minimums[element] = Int(count)
        end
        intervals = Dict{Tuple{String,String},Tuple{Float64,Float64}}()
        for (pair,interval) in ratios
            pair isa Tuple && length(pair)==2 || throw(ArgumentError("Ratio keys must be (numerator,denominator) tuples."))
            key = (String(pair[1]),String(pair[2]))
            all(e->haskey(ELEMENTS,e) && e!="e",key) || throw(ArgumentError("Unknown ratio element."))
            lo,hi = checked_interval(interval,"ratio")
            lo >= 0 || throw(ArgumentError("Ratio limits must be nonnegative."))
            intervals[key] = (lo,hi)
        end
        dbe_interval = isnothing(dbe) ? nothing : checked_interval(dbe,"DBE")
        new(minimums,intervals,dbe_interval,valence)
    end
end

function checked_interval(interval,label)
    (interval isa Tuple || interval isa AbstractVector) && length(interval)==2 || throw(ArgumentError("$label must be a (minimum,maximum) interval."))
    lo,hi = interval
    lo isa Real && hi isa Real && !isnan(lo) && !isnan(hi) && lo<=hi || throw(ArgumentError("Invalid $label interval."))
    return Float64(lo),Float64(hi)
end

"Counts after atom-pool and atom-loss availability bounds, with first rejection reasons."
mutable struct SearchDiagnostics
    mass_matches::Int
    rejections::Dict{Symbol,Int}
    status::Symbol
end
SearchDiagnostics() = SearchDiagnostics(0,Dict{Symbol,Int}(),:no_mass_match)
reject!(diagnostics,reason) = (diagnostics.rejections[reason] = get(diagnostics.rejections,reason,0)+1; nothing)

"Matching formulas, diagnostic counts, and a snapshot of the search parameters."
struct SearchResult
    matches::Vector{Compound}
    diagnostics::SearchDiagnostics
    settings::NamedTuple
end
