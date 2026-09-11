# JSON-safe snapshots. Infinite resolving power and open constraint bounds are
# strings, so snapshots remain valid JSON rather than nonstandard numeric tokens.
portable_value(x::Union{Nothing,Bool,Integer,AbstractString}) = x
portable_value(x::Symbol) = String(x)
portable_value(x::AbstractFloat) = isfinite(x) ? x : string(x)
portable_value(x::Real) = portable_value(Float64(x))
portable_value(x::Union{Tuple,AbstractVector}) = [portable_value(v) for v in x]
portable_value(x::NamedTuple) = Dict(String(k)=>portable_value(v) for (k,v) in pairs(x))
portable_value(x::AbstractDict) = Dict(String(k)=>portable_value(v) for (k,v) in x)
function portable_value(x::SearchConstraints)
    return Dict("min_atoms"=>portable_value(x.min_atoms),"valence"=>String(x.valence),"dbe"=>portable_value(x.dbe),
                "ratios"=>[Dict("numerator"=>a,"denominator"=>b,"minimum"=>portable_value(lo),"maximum"=>portable_value(hi))
                            for ((a,b),(lo,hi)) in sort!(collect(x.ratios);by=first)])
end
function portable_value(x::Union{Compound,IonSpec,IsotopePopulation,IsotopeContribution,IsotopePeak,PatternResult,
                                  SearchDiagnostics,SearchResult,PeakMatch,RankedCandidate})
    return Dict(String(field)=>portable_value(getfield(x,field)) for field in fieldnames(typeof(x)))
end

"""
    result_dict(result)

Return a JSON-safe schema-version-2 snapshot of a PatternResult, SearchResult or
RankedCandidate. Includes inputs, settings and isotope-data fingerprint. Infinite
bounds/resolving power are represented as "Inf" or "-Inf". Serialize with JSON.jl;
this is an export format, not an automatic result-restoration API.
"""
function result_dict(result::Union{PatternResult,SearchResult,RankedCandidate})
    kind = result isa PatternResult ? "pattern" : result isa SearchResult ? "search" : "ranked_candidate"
    return Dict("schema_version"=>2,"kind"=>kind,"data"=>portable_value(result))
end
