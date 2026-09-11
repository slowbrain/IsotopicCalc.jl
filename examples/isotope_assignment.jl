# Run: julia --project=. examples/isotope_assignment.jl [output.json]
using IsotopicCalc
using JSON

# Synthetic centroids for a reproducible example, with a small mass offset.
reference = isotopic_pattern("C2H4O";R=1000,print=false)
observed = [(mass+0.00003,1000*abundance) for (mass,abundance) in reference]
constraints = SearchConstraints(min_atoms=Dict("O"=>1),dbe=(0,5))
ranked = rank_formulas(first(observed)[1],observed;
    search_options=(atom_pool=Dict("C"=>3,"H"=>8,"O"=>3),tolerance_ppm=1000,constraints=constraints),
    pattern_options=(R=1000,),match_tolerance_da=0.05)
for candidate in ranked
    println(candidate.compound.formula, "  ppm=",round(candidate.compound.ppm;digits=2),
            "  isotope similarity=",round(candidate.isotope_similarity;digits=6),
            "  missing=",length(candidate.missing_peaks),"  unexplained=",length(candidate.unexplained_peaks))
end
if !isempty(ARGS) && !isempty(ranked)
    open(ARGS[1],"w") do io
        JSON.print(io,result_dict(first(ranked)),2)
    end
end
