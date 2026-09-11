"One measured/predicted peak pair, with signed observed-minus-predicted residuals."
struct PeakMatch
    observed_index::Int
    predicted_index::Int
    mass_error_da::Float64
    mass_error_ppm::Float64
    observed_intensity::Float64
    fitted_intensity::Float64
    intensity_residual::Float64
end

"Candidate evidence; scores are descriptive similarities, not identification probabilities."
struct RankedCandidate
    compound::Compound
    mass_error_da::Float64
    mass_score::Float64
    isotope_similarity::Float64
    intensity_scale::Float64
    matches::Vector{PeakMatch}
    missing_peaks::Vector{Int}
    unexplained_peaks::Vector{Int}
    pattern::PatternResult
    match_settings::NamedTuple
end

function observed_peaks(observed)
    values = Tuple{Float64,Float64}[]
    for peak in observed
        (peak isa Tuple || peak isa AbstractVector) && length(peak)==2 || throw(ArgumentError("Observed peaks must be (mass,intensity) pairs."))
        mass,intensity = peak
        mass isa Real && isfinite(mass) && mass>0 || throw(ArgumentError("Observed masses must be finite and positive."))
        intensity isa Real && isfinite(intensity) && intensity>=0 || throw(ArgumentError("Observed intensities must be finite and nonnegative."))
        push!(values,(Float64(mass),Float64(intensity)))
    end
    isempty(values) && throw(ArgumentError("Observed isotope cluster cannot be empty."))
    maximum(last,values)>0 || throw(ArgumentError("Observed intensities cannot all be zero."))
    length(unique(first.(values))) == length(values) || throw(ArgumentError("Observed centroid masses must be unique."))
    return values
end

# Ordered one-to-one alignment: maximize the number of admissible pairs, then
# minimize their summed squared mass error normalized by the matching window.
function align_peaks(predicted,observed,ppm,da)
    order = sortperm(observed;by=first)
    n,m = length(predicted),length(order)
    counts = zeros(Int,n+1,m+1)
    errors = zeros(Float64,n+1,m+1)
    steps = zeros(UInt8,n+1,m+1)
    for i in 1:n, j in 1:m
        count,cost,step = counts[i,j+1],errors[i,j+1],UInt8(1)
        if counts[i+1,j]>count || (counts[i+1,j]==count && errors[i+1,j]<cost)
            count,cost,step = counts[i+1,j],errors[i+1,j],UInt8(2)
        end
        mass = predicted[i].mass
        difference = observed[order[j]][1]-mass
        window = max(da,mass*ppm/1e6)
        if abs(difference)<=window
            candidate_count = counts[i,j]+1
            candidate_error = errors[i,j]+(window==0 ? 0.0 : (difference/window)^2)
            if candidate_count>count || (candidate_count==count && candidate_error<=cost)
                count,cost,step = candidate_count,candidate_error,UInt8(3)
            end
        end
        counts[i+1,j+1],errors[i+1,j+1],steps[i+1,j+1] = count,cost,step
    end
    pairs = Tuple{Int,Int}[]
    i,j = n,m
    while i>0 && j>0
        step = steps[i+1,j+1]
        if step==3
            push!(pairs,(i,order[j])); i-=1; j-=1
        elseif step==1
            i-=1
        else
            j-=1
        end
    end
    reverse!(pairs)
    return pairs
end

"""
    rank_formulas(mz_input, observed; search_options=(;), pattern_options=(;),
                  match_tolerance_ppm=10, match_tolerance_da=0.001)

Search formulas and compare their predicted isotope clusters with measured
`(mass,intensity)` centroids. Supply one local cluster, including missing/extra
peak evidence, consistently for all candidates. Peaks outside the supplied
cluster cannot provide evidence. Input indices in results refer to the original
observed order. The search's ion/mode determines each prediction.

One-to-one mass matching maximizes matched peak count then minimizes normalized
mass error. Cosine intensity similarity includes all predicted and observed
peaks (unmatched entries contribute zero to the dot product). A least-squares
intensity scale uses matched observations, treating missing predictions as zero.
Sort by descending isotope similarity, then absolute precursor ppm error.
`mass_score = exp(-0.5*(mass_error/max(absolute_tolerance,calculated_mz*ppm/1e6))^2)`;
neither score is a calibrated identification probability.
"""
function rank_formulas(mz_input::Real,observed;search_options::NamedTuple=(;),pattern_options::NamedTuple=(;),
                       match_tolerance_ppm::Real=10,match_tolerance_da::Real=0.001)
    isfinite(match_tolerance_ppm) && match_tolerance_ppm>=0 || throw(ArgumentError("Matching ppm tolerance must be finite and nonnegative."))
    isfinite(match_tolerance_da) && match_tolerance_da>=0 || throw(ArgumentError("Matching absolute tolerance must be finite and nonnegative."))
    values = observed_peaks(observed)
    for key in (:ion,:annotated,:print)
        haskey(pattern_options,key) && throw(ArgumentError("'$key' is controlled by rank_formulas; configure the ion through search_options."))
    end
    search = find_formula(mz_input;merge(search_options,(print=false,diagnostics=true))...)
    maximum_intensity = maximum(last,values)
    intensities = [p[2]/maximum_intensity for p in values]
    observed_norm = sqrt(sum(abs2,intensities))
    results = RankedCandidate[]
    for compound in search.matches
        pattern = isotopic_pattern(compound.formula;merge(pattern_options,
                    (ion=search.settings.ion,annotated=true,print=false))...)
        pairs = align_peaks(pattern.peaks,values,match_tolerance_ppm,match_tolerance_da)
        predicted = [p.abundance for p in pattern.peaks]
        dot_product = sum(Float64[predicted[i]*intensities[j] for (i,j) in pairs])
        predicted_squared = sum(abs2,predicted)
        similarity = clamp(dot_product/(sqrt(predicted_squared)*observed_norm),0.0,1.0)
        fitted_scale = dot_product/predicted_squared*maximum_intensity
        matches = PeakMatch[]
        for (i,j) in pairs
            error = values[j][1]-pattern.peaks[i].mass
            fitted = predicted[i]*fitted_scale
            push!(matches,PeakMatch(j,i,error,error/pattern.peaks[i].mass*1e6,values[j][2],fitted,values[j][2]-fitted))
        end
        matched_predicted = Set(i for (i,j) in pairs)
        matched_observed = Set(j for (i,j) in pairs)
        missing = [i for i in eachindex(predicted) if !(i in matched_predicted)]
        unexplained = [j for j in eachindex(values) if !(j in matched_observed)]
        error = Float64(mz_input-compound.mz)
        window = max(search.settings.tolerance_da,compound.mz*search.settings.tolerance_ppm/1e6)
        mass_score = window==0 ? (error==0 ? 1.0 : 0.0) : exp(-0.5*(error/window)^2)
        settings = (match_tolerance_ppm=match_tolerance_ppm,match_tolerance_da=match_tolerance_da,
                    search=search.settings,observed=copy(values))
        push!(results,RankedCandidate(compound,error,mass_score,similarity,fitted_scale,matches,missing,unexplained,pattern,settings))
    end
    sort!(results;by=result->(-result.isotope_similarity,abs(result.compound.ppm),result.compound.formula))
    return results
end
