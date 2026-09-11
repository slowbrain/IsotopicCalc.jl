"An isotopologue contributing absolute probability to a resolved/merged peak."
struct IsotopeContribution
    composition::Dict{String,Int}
    mass::Float64
    probability::Float64
    population::Int
end

"A centroid, its base-peak-relative abundance, and absolute probability."
struct IsotopePeak
    mass::Float64
    abundance::Float64
    probability::Float64
    contributors::Vector{IsotopeContribution}
end

"""
    IsotopePopulation(formula, weight; abundances=Dict())

A molecular population in a mixture. All population weights must sum to one;
all formulas must have the same elemental counts as the base formula.
"""
struct IsotopePopulation
    formula::String
    weight::Float64
    abundances::Dict{String,Dict{String,Float64}}
    function IsotopePopulation(formula::AbstractString,weight::Real;abundances=Dict())
        isfinite(weight) && 0 <= weight <= 1 || throw(ArgumentError("Population weight must be in [0,1]."))
        for atom in keys(parse_formula(formula)); principal_mass(atom); end
        new(String(formula),Float64(weight),validate_abundances(abundances))
    end
end

"Annotated pattern, probability accounting and a snapshot of calculation settings."
struct PatternResult
    formula::String
    ion_composition::Dict{String,Int}
    charge::Int
    peaks::Vector{IsotopePeak}
    retained_probability::Float64
    cutoff_probability::Float64
    pruned_probability::Float64
    pruning_bound::Float64
    settings::NamedTuple
end

function validate_abundances(overrides)
    overrides isa AbstractDict || throw(ArgumentError("abundances must map elements to isotope-probability dictionaries."))
    result = Dict{String,Dict{String,Float64}}()
    for (key,values) in overrides
        element = String(key)
        haskey(ELEMENTS,element) && !(element in ("e","D")) || throw(ArgumentError("Unknown or fixed-isotope abundance element '$element'."))
        values isa AbstractDict || throw(ArgumentError("Abundances for '$element' must be a dictionary."))
        probabilities = Dict{String,Float64}()
        for (isotope,p) in values
            label = String(isotope)
            label in ELEMENTS[element]["Symbol"] || throw(ArgumentError("Unknown isotope '$label' for '$element'."))
            p isa Real && isfinite(p) && 0 <= p <= 1 || throw(ArgumentError("Isotope probabilities must be finite and in [0,1]."))
            probabilities[label] = Float64(p)
        end
        total = sum(Float64[p for p in Base.values(probabilities)])
        isapprox(total,1.0;atol=1e-12,rtol=0) || throw(ArgumentError("Isotope probabilities for '$element' must sum to one; omitted isotopes have zero probability."))
        for isotope in keys(probabilities); probabilities[isotope] /= total; end
        result[element] = probabilities
    end
    return result
end

function population_isotopes(atom, overrides)
    masses, probabilities = isotope_data(atom)
    labels = haskey(ELEMENTS,atom) ? String.(ELEMENTS[atom]["Symbol"]) : [atom]
    atom == "D" && (labels = ["2H"])
    if haskey(overrides,atom)
        probabilities = [get(overrides[atom],label,0.0) for label in labels]
    end
    return masses,probabilities,labels
end

function elemental_counts(formula)
    result = Dict{String,Int}()
    for (atom,count) in parse_formula(formula)
        base = atom == "D" ? "H" : haskey(ELEMENTS,atom) ? atom : get_base_element(atom)
        add_atom_count!(result,base,count)
    end
    return result
end

function pattern_populations(formula,overrides,mixture)
    isnothing(mixture) && return [IsotopePopulation(formula,1.0;abundances=overrides)]
    isempty(overrides) || throw(ArgumentError("Use per-population abundances when supplying a mixture."))
    populations = IsotopePopulation[]
    for population in mixture
        population isa IsotopePopulation || throw(ArgumentError("mixture must contain IsotopePopulation values."))
        push!(populations,IsotopePopulation(population.formula,population.weight;abundances=population.abundances))
    end
    isempty(populations) && throw(ArgumentError("Mixture cannot be empty."))
    isapprox(sum(p.weight for p in populations),1.0;atol=1e-12,rtol=0) || throw(ArgumentError("Population weights must sum to one."))
    expected = elemental_counts(formula)
    all(p -> elemental_counts(p.formula) == expected,populations) || throw(ArgumentError("Population formulas must have the same elemental counts as the base formula."))
    total = sum(p.weight for p in populations)
    return [IsotopePopulation(p.formula,p.weight/total;abundances=p.abundances) for p in populations]
end

function prune_contributions(contributions,budget)
    budget <= 0 && return contributions
    retained = trues(length(contributions))
    discarded = 0.0
    for index in sortperm(contributions;by=c->c.probability)
        discarded + contributions[index].probability > budget && break
        discarded += contributions[index].probability
        retained[index] = false
    end
    return contributions[retained]
end

function annotated_element(masses,probabilities,labels,count,budget,::Val{N}) where N
    states = isotope_states(probabilities,count,budget,Val(N))
    return [IsotopeContribution(Dict(labels[i]=>composition[i] for i in 1:N if composition[i]>0),
             sum(composition[i]*masses[i] for i in 1:N),probability,1) for (composition,probability) in states]
end

function annotated_convolve(left,right,budget)
    result = IsotopeContribution[]
    for a in left, b in right
        a.probability*b.probability == 0 && continue
        composition = copy(a.composition)
        for (atom,count) in b.composition; add_atom_count!(composition,atom,count); end
        push!(result,IsotopeContribution(composition,a.mass+b.mass,a.probability*b.probability,1))
    end
    return prune_contributions(result,budget)
end

function annotated_groups(contributions,R)
    groups = Vector{IsotopeContribution}[]
    anchor = -Inf
    for contribution in sort(contributions;by=c->c.mass)
        if isempty(groups) || contribution.mass-anchor > anchor/R
            push!(groups,IsotopeContribution[])
            anchor = contribution.mass
        end
        push!(last(groups),contribution)
    end
    return groups
end

function calculate_pattern(formula,spec,cutoff,R,overrides,mixture,annotated)
    populations = pattern_populations(formula,overrides,mixture)
    error_budget = cutoff == 0 ? 0.0 : min(1e-12,cutoff*1e-8)
    tuples = Tuple{Float64,Float64}[]
    contributions = IsotopeContribution[]
    for (population_id,population) in enumerate(populations)
        composition = ion_composition(population.formula,spec)
        population.weight == 0 && continue
        atoms = sort!(collect(keys(composition)))
        steps = sum(Float64[composition[atom] for atom in atoms])+length(atoms)
        budget = error_budget/max(steps,1)
        distribution = [(0.0,1.0)]
        details = [IsotopeContribution(Dict{String,Int}(),0.0,1.0,1)]
        for atom in atoms
            masses,probabilities,labels = population_isotopes(atom,population.abundances)
            if annotated
                details = annotated_convolve(details,annotated_element(masses,probabilities,labels,composition[atom],budget,Val(length(masses))),budget)
            else
                distribution = convolve(distribution,element_distribution(masses,probabilities,composition[atom],budget,Val(length(masses))),budget)
            end
        end
        electron_adjustment = spec.charge*principal_mass("e")
        divisor = max(abs(spec.charge),1)
        if annotated
            for c in details
                push!(contributions,IsotopeContribution(c.composition,(c.mass-electron_adjustment)/divisor,c.probability*population.weight,population_id))
            end
        else
            append!(tuples,[((m-electron_adjustment)/divisor,p*population.weight) for (m,p) in distribution])
        end
    end
    if !annotated
        all(p->p[1]>0,tuples) || throw(ArgumentError("Ion masses must be positive."))
        grouped = normalize_abundances(group_by_resolution(tuples,R))
        return filter(p->p[2]>=cutoff,grouped)
    end
    all(c->c.mass>0,contributions) || throw(ArgumentError("Ion masses must be positive."))
    groups = annotated_groups(contributions,R)
    probabilities = [sum(c.probability for c in group) for group in groups]
    maximum_probability = maximum(probabilities)
    peaks = IsotopePeak[]
    removed = 0.0
    for (group,probability) in zip(groups,probabilities)
        abundance = probability/maximum_probability
        if abundance >= cutoff
            mass = sum(c.mass*c.probability for c in group)/probability
            push!(peaks,IsotopePeak(mass,abundance,probability,group))
        else
            removed += probability
        end
    end
    settings = (R=R,abundance_cutoff=cutoff,ion=deepcopy(spec),populations=deepcopy(populations),
                package_version=string(parentmodule(@__MODULE__).PACKAGE_VERSION),
                isotope_source="NIST Atomic Weights and Isotopic Compositions (bundled snapshot)",
                isotope_data_sha256=parentmodule(@__MODULE__).ISOTOPE_DATA_SHA256)
    retained = sum(Float64[p.probability for p in peaks])
    return PatternResult(sum_formula(parse_formula(formula)),ion_composition(formula,spec),spec.charge,peaks,
                         retained,removed,max(0.0,1-sum(probabilities)),error_budget,settings)
end
