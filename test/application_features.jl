@testset "Explicit ions" begin
    h,o,n,cl,e = 1.00782503223,15.99491461957,14.00307400443,34.968852682,0.0005485
    acetone = 36+6h+o
    for (notation,expected) in (("[M+NH4]+",acetone+n+4h-e),
                               ("[M+Cl]-",acetone+cl+e),("[M+H-H2O]+",acetone-h-o-e),
                               ("[M+2H]2+",(acetone+2h-2e)/2),("[M-H]-",acetone-h+e),("[M]2+",(acetone-2e)/2))
        @test ion_mass("C3H6O";ion=notation) ≈ expected atol=1e-11
        @test isotopic_pattern("C3H6O";ion=notation,print=false)[1][1] ≈ expected atol=1e-11
    end
    @test ion_mass("C3H6O";ion=IonSpec(add="H",remove="H2O",charge=1)) ≈ acetone-h-o-e
    @test ion_mass("C3H7O";ion=IonSpec(charge=1)) ≈ acetone+h-e
    @test ion_mass("Fe";ion="[M]+") ≈ 55.93493633-e
    @test ion_mass("C3H6O";ion="[M+[13C]]+") ≈ acetone+13.00335483507-e
    @test_throws ArgumentError ion_mass("CO2";ion="[M-H2O]-")
    @test_throws ArgumentError ion_mass("H2O";ion=nothing)
    for bad in ("[M+H]0+","[M+] +","[M+H-H2O]2","[M++H]+","[M+0H]+","Xx+","[M+Xx]+")
        @test_throws ArgumentError IonSpec(bad)
    end
    original = Dict("H"=>1)
    spec = IonSpec(original,Dict{String,Int}(),1)
    original["H"] = 2
    @test spec.additions["H"] == 1
end

@testset "Annotated patterns, populations and enrichment" begin
    enriched = Dict("C"=>Dict("12C"=>0.8,"13C"=>0.2))
    result = isotopic_pattern("C2";abundances=enriched,annotated=true,R=Inf,print=false)
    @test result isa PatternResult
    @test [p.probability for p in result.peaks] ≈ [0.64,0.32,0.04]
    @test [p.abundance for p in result.peaks] ≈ [1,0.5,0.0625]
    @test result.peaks[2].contributors[1].composition == Dict("12C"=>1,"13C"=>1)
    @test result.retained_probability ≈ 1
    @test result.cutoff_probability ≈ 0 atol=1e-12
    @test result.pruned_probability <= result.pruning_bound+1e-12
    @test length(result.settings.isotope_data_sha256)==64
    @test result.settings.package_version == "0.8.0"
    enriched["C"]["13C"] = 0.9
    @test result.settings.populations[1].abundances["C"]["13C"] == 0.2
    @test IsotopicCalc.ELEMENTS["C"]["Isotopic Composition"][2] == 0.0107
    overrides = Dict("C"=>Dict("12C"=>0.8,"13C"=>0.2))
    full_label = isotopic_pattern("[13C]2";abundances=overrides,annotated=true,print=false)
    @test length(full_label.peaks)==1
    @test full_label.peaks[1].probability==1
    @test full_label.peaks[1].mass ≈ 26.00670967014
    cut = isotopic_pattern("C2";abundances=overrides,abundance_cutoff=0.1,annotated=true,print=false)
    @test cut.retained_probability ≈ 0.96
    @test cut.cutoff_probability ≈ 0.04
    mixture = [IsotopePopulation("C2",0.75;abundances=Dict("C"=>Dict("12C"=>1.0))),IsotopePopulation("[13C]2",0.25)]
    mixed = isotopic_pattern("C2";mixture=mixture,annotated=true,print=false)
    @test [p.probability for p in mixed.peaks] ≈ [0.75,0.25]
    @test [only(p.contributors).population for p in mixed.peaks] == [1,2]
    @test [p.mass for p in mixed.peaks] ≈ [24,26.00670967014]
    @test isempty(filter(p->25<p.mass<26,mixed.peaks)) # Mixtures are not independent enrichment.
    for formula in ("C3H6O","C100","[13C]C2D6O","Fe","H2O")
        for ion in ("[M]","[M+H]+","[M+Cl]-")
            detailed = isotopic_pattern(formula;ion=ion,annotated=true,print=false)
            plain = isotopic_pattern(formula;ion=ion,print=false)
            @test length(plain)==length(detailed.peaks)
            @test first.(plain) ≈ [p.mass for p in detailed.peaks] atol=1e-10
            @test last.(plain) ≈ [p.abundance for p in detailed.peaks] atol=1e-10
            @test detailed.retained_probability+detailed.cutoff_probability+detailed.pruned_probability ≈ 1 atol=1e-12
            for peak in detailed.peaks
                @test sum(c.probability for c in peak.contributors) ≈ peak.probability
                @test sum(c.mass*c.probability for c in peak.contributors)/peak.probability ≈ peak.mass
            end
        end
    end
    for bad in (Dict("C"=>Dict("13C"=>0.2)),Dict("C"=>Dict("14C"=>1.0)),Dict("C"=>Dict("13C"=>NaN)),Dict("D"=>Dict("2D"=>1.0)),Dict("C"=>[0.8,0.2]))
        @test_throws ArgumentError isotopic_pattern("C2";abundances=bad,print=false)
    end
    @test_throws ArgumentError IsotopePopulation("C2",-0.1)
    @test_throws ArgumentError isotopic_pattern("C2";mixture=[IsotopePopulation("C",1)],print=false)
    @test_throws ArgumentError isotopic_pattern("C2";mixture=[IsotopePopulation("C2",0.5)],print=false)
    @test_throws ArgumentError isotopic_pattern("C2";mixture=mixture,abundances=overrides,print=false)
end

@testset "Search constraints, diagnostics and ion compositions" begin
    acetone = 58.04186481295
    pool = Dict("C"=>3,"H"=>8,"O"=>2)
    report = find_formula(acetone;atom_pool=pool,tolerance_ppm=1,diagnostics=true,print=false)
    @test report isa SearchResult
    @test report.diagnostics.status == :matched
    @test report.diagnostics.mass_matches == length(report.matches)
    @test only(report.matches).formula == "C3H6O"
    constraints = SearchConstraints(min_atoms=Dict("O"=>2))
    filtered = find_formula(acetone;atom_pool=pool,tolerance_ppm=1,constraints=constraints,diagnostics=true,print=false)
    @test isempty(filtered.matches)
    @test filtered.diagnostics.status == :filtered_out
    @test filtered.diagnostics.rejections[:minimum_atoms]==1
    constraints.min_atoms["O"]=1
    @test filtered.settings.constraints.min_atoms["O"]==2
    no_mass = find_formula(acetone+0.01;atom_pool=pool,tolerance_ppm=1,diagnostics=true,print=false)
    @test no_mass.diagnostics.status == :no_mass_match
    @test no_mass.diagnostics.mass_matches==0
    absolute = find_formula(acetone+0.001;atom_pool=pool,tolerance_ppm=0,tolerance_da=0.002,print=false)
    @test any(c->c.formula=="C3H6O",absolute)
    @test isempty(find_formula(acetone+0.003;atom_pool=pool,tolerance_ppm=0,tolerance_da=0.002,print=false))
    for dbe in ((0,0),(2,4))
        r = find_formula(acetone;atom_pool=pool,tolerance_ppm=1,constraints=SearchConstraints(dbe=dbe),diagnostics=true,print=false)
        @test r.diagnostics.rejections[:dbe]==1
    end
    @test !isempty(find_formula(acetone;atom_pool=pool,tolerance_ppm=1,constraints=SearchConstraints(dbe=(1,1)),print=false))
    cs2 = 75.9441423488
    r = find_formula(cs2;atom_pool=Dict("C"=>1,"S"=>2),tolerance_ppm=1,
                     constraints=SearchConstraints(ratios=Dict(("S","C")=>(0,2))),print=false)
    @test only(r).formula=="CS2"
    r = find_formula(acetone;atom_pool=pool,tolerance_ppm=1,
                    constraints=SearchConstraints(ratios=Dict(("H","C")=>(0,1))),diagnostics=true,print=false)
    @test r.diagnostics.rejections[:ratios]==1
    protonated = acetone+1.00782503223-0.0005485
    direct = find_formula(protonated;atom_pool=pool,tolerance_ppm=1,mode=:ion,ion=IonSpec(charge=1),print=false)
    @test only(direct).formula=="C3H7O"
    acetate = 24+3*1.00782503223+2*15.99491461957+0.0005485
    direct = find_formula(acetate;atom_pool=pool,tolerance_ppm=1,mode=:ion,ion=IonSpec(charge=-1),print=false)
    @test only(direct).formula=="C2H3O2"
    radical = find_formula(acetone-0.0005485;atom_pool=pool,tolerance_ppm=1,mode=:ion,
                          ion=IonSpec(charge=1),constraints=SearchConstraints(valence=:none),print=false)
    @test only(radical).formula=="C3H6O"
    for notation in ("[M+NH4]+","[M+Cl]-","[M+H-H2O]+")
        r = find_formula(ion_mass("C3H6O";ion=notation);atom_pool=pool,ion=notation,tolerance_ppm=1,print=false)
        @test any(c->c.formula=="C3H6O",r)
    end
    infeasible = find_formula(40;ion="[M-H]-",atom_pool=Dict("C"=>4),diagnostics=true,print=false)
    @test infeasible.diagnostics.status==:infeasible_ion
    @test_throws ArgumentError find_formula(protonated;mode=:ion,ion="[M+H]+",print=false)
    @test_throws ArgumentError find_formula(10;tolerance_da=-1)
    @test_throws ArgumentError SearchConstraints(min_atoms=Dict("H"=>-1))
    @test_throws ArgumentError SearchConstraints(ratios=Dict(("H","C")=>(2,1)))
    @test_throws ArgumentError SearchConstraints(dbe=(NaN,2))
end

# Independent multinomial reference grouped by nominal isotope increment.
function observed_acetaldehyde()
    groups = Dict{Int,Tuple{Float64,Float64}}()
    for c in 0:2, h in 0:4, (increment,mass_o,p_o) in ((0,15.99491461957,0.99757),(1,16.99913175650,0.00038),(2,17.99915961286,0.00205))
        probability = binomial(2,c)*0.0107^c*0.9893^(2-c)*binomial(4,h)*0.000115^h*0.999885^(4-h)*p_o
        mass = 24+c*1.00335483507+(4-h)*1.00782503223+h*2.01410177812+mass_o
        weighted,total = get(groups,c+h+increment,(0.0,0.0))
        groups[c+h+increment]=(weighted+mass*probability,total+probability)
    end
    maximum_probability = maximum(last,values(groups))
    return sort([(weighted/probability,1000*probability/maximum_probability) for (weighted,probability) in values(groups) if probability/maximum_probability>=1e-5])
end

@testset "Observed-pattern ranking and evidence" begin
    observed = observed_acetaldehyde()
    search = (atom_pool=Dict("C"=>2,"H"=>4,"O"=>2),tolerance_ppm=1000)
    ranked = rank_formulas(43.98982923914,observed;search_options=search,pattern_options=(R=1000,),match_tolerance_da=0.05)
    @test first(ranked).compound.formula=="C2H4O"
    @test first(ranked).isotope_similarity ≈ 1 atol=1e-10
    @test first(ranked).intensity_scale ≈ 1000 atol=1e-6
    @test isempty(first(ranked).missing_peaks)
    @test isempty(first(ranked).unexplained_peaks)
    co2 = only(filter(c->c.compound.formula=="CO2",ranked))
    @test co2.mass_score>first(ranked).mass_score
    @test co2.isotope_similarity<first(ranked).isotope_similarity
    @test maximum(abs(m.intensity_residual) for m in first(ranked).matches)<1e-6
    @test maximum(abs(m.mass_error_da) for m in first(ranked).matches)<1e-9
    reversed = rank_formulas(43.98982923914,reverse(observed);search_options=search,pattern_options=(R=1000,),match_tolerance_da=0.05)
    @test first(reversed).isotope_similarity ≈ first(ranked).isotope_similarity
    @test first(reversed).matches[1].observed_index==length(observed)
    incomplete = [observed[1],(observed[1][1]+0.1,100.0)]
    damaged = first(rank_formulas(44.02621474849,incomplete;search_options=(atom_pool=Dict("C"=>2,"H"=>4,"O"=>1),tolerance_ppm=1),pattern_options=(R=1000,),match_tolerance_da=0.005))
    @test damaged.isotope_similarity<1
    @test !isempty(damaged.missing_peaks)
    @test damaged.unexplained_peaks==[2]
    @test length(unique(m.observed_index for m in damaged.matches))==length(damaged.matches)
    @test length(unique(m.predicted_index for m in damaged.matches))==length(damaged.matches)
    # Charge-only actual-ion search must use that composition once, without an extra proton.
    mz = 36+7*1.00782503223+15.99491461957-0.0005485
    ions = rank_formulas(mz,[(mz,100.0)];search_options=(atom_pool=Dict("C"=>3,"H"=>7,"O"=>1),tolerance_ppm=1,mode=:ion,ion=IonSpec(charge=1)))
    @test only(ions).compound.formula=="C3H7O"
    @test only(ions).pattern.peaks[1].mass ≈ mz
    for bad in ([],[(1.0,0.0)],[(1.0,-1.0)],[(NaN,1.0)],[(1.0,1.0),(1.0,2.0)])
        @test_throws ArgumentError rank_formulas(12,bad)
    end
    @test_throws ArgumentError rank_formulas(12,[(12,1)];pattern_options=(ion=IonSpec(),))
    @test_throws ArgumentError rank_formulas(12,[(12,1)];match_tolerance_da=-1)
end

@testset "Portable result snapshots" begin
    pattern = isotopic_pattern("C2";R=Inf,annotated=true,print=false)
    snapshot = result_dict(pattern)
    @test snapshot["schema_version"]==2
    @test snapshot["kind"]=="pattern"
    @test snapshot["data"]["settings"]["R"]=="Inf"
    decoded = IsotopicCalc.JSON.parse(IsotopicCalc.JSON.json(snapshot))
    @test decoded["data"]["peaks"][1]["contributors"][1]["composition"]["12C"]==2
    report = find_formula(24;diagnostics=true,print=false,constraints=SearchConstraints(ratios=Dict(("H","C")=>(0,Inf))))
    decoded = IsotopicCalc.JSON.parse(IsotopicCalc.JSON.json(result_dict(report)))
    @test decoded["data"]["settings"]["constraints"]["ratios"][1]["maximum"]=="Inf"
    rank = only(rank_formulas(24,[(24.0,100.0)];search_options=(atom_pool=Dict("C"=>2),tolerance_ppm=1)))
    decoded = IsotopicCalc.JSON.parse(IsotopicCalc.JSON.json(result_dict(rank)))
    @test decoded["data"]["compound"]["formula"]=="C2"
    @test decoded["data"]["match_settings"]["observed"][1][2]==100
end

@testset "Absolute/ppm windows and constraints against enumeration" begin
    pool = Dict("C"=>2,"H"=>6,"O"=>2)
    constraints = SearchConstraints(min_atoms=Dict("H"=>2),dbe=(0,2),ratios=Dict(("H","C")=>(0,2)))
    for (charge,ion,delta) in ((0,"[M]",0.0),(2,"[M+2H]2+",2*1.00782503223),(-1,"[M-H]-",-1.00782503223))
        for (ppm,da) in ((0.0,0.01),(100.0,0.0),(10000.0,0.005),(100.0,2.0))
            target = (28.0312+delta-charge*0.0005485)/max(abs(charge),1)
            expected = Set{NTuple{3,Int}}()
            for c in 0:2,h in 0:6,o in 0:2
                charge<0 && h==0 && continue
                h>=2 || continue
                dbe = (2c+2-h)/2
                isinteger(dbe) && 0<=dbe<=2 || continue
                c>0 && h/c>2 && continue
                mz = (12c+1.00782503223h+15.99491461957o+delta-charge*0.0005485)/max(abs(charge),1)
                (abs((target-mz)/mz*1e6)<=ppm || abs(target-mz)<=da) && push!(expected,(c,h,o))
            end
            r = find_formula(target;atom_pool=pool,ion=ion,tolerance_ppm=ppm,tolerance_da=da,
                             constraints=constraints,diagnostics=true,print=false,heuristics=false)
            actual = Set(Tuple(get(IsotopicCalc.IsotopicPatternModule.parse_formula(x.formula),e,0) for e in ("C","H","O")) for x in r.matches)
            @test actual==expected
            @test length(r.matches)+sum(Int[v for v in values(r.diagnostics.rejections)])==r.diagnostics.mass_matches
        end
    end
end
