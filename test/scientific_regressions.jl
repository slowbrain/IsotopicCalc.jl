# Independent reference values and binomial/multinomial probabilities. These
# tests must not derive expected masses or abundances through the package API.
@testset "Scientific correctness regressions" begin
    @testset "Carbon isotope probabilities" begin
        p = 0.0107
        for n in (5, 10, 100)
            probabilities = [Float64(binomial(big(n), k)) * p^k * (1-p)^(n-k) for k in 0:n]
            relative = probabilities ./ maximum(probabilities)
            expected = [(12.0*n + k*1.00335483507, relative[k+1]) for k in 0:n if relative[k+1] >= 1e-5]
            actual = isotopic_pattern("C$n"; R=Inf, print=false)
            @test length(actual) == length(expected)
            for (observed, reference) in zip(actual, expected)
                @test observed[1] ≈ reference[1] atol=1e-10 rtol=0
                @test observed[2] ≈ reference[2] atol=1e-11 rtol=0
            end
        end
        c2 = isotopic_pattern("C2"; abundance_cutoff=0.015, R=Inf, print=false)
        @test length(c2) == 2
        @test c2[2][2] ≈ 2p/(1-p) atol=1e-14
        cl2 = isotopic_pattern("Cl2"; abundance_cutoff=0.1, R=Inf, print=false)
        @test length(cl2) == 3
        @test cl2[3][2] ≈ (0.2424/0.7576)^2 atol=1e-14
    end

    @testset "Mixed element probabilities" begin
        # Acetone: independently enumerate counts of 13C and 2H, and the
        # three possible oxygen isotopes, using closed-form probabilities.
        reference = Tuple{Float64,Float64}[]
        for c in 0:3, h in 0:6, (o_mass, o_probability) in
                ((15.99491461957, 0.99757), (16.99913175650, 0.00038), (17.99915961286, 0.00205))
            probability = binomial(3,c)*0.0107^c*0.9893^(3-c) *
                          binomial(6,h)*0.000115^h*0.999885^(6-h) * o_probability
            mass = 36 + c*1.00335483507 + (6-h)*1.00782503223 + h*2.01410177812 + o_mass
            push!(reference, (mass,probability))
        end
        max_probability = maximum(last, reference)
        expected = sort([(m,p/max_probability) for (m,p) in reference if p/max_probability >= 1e-5])
        actual = isotopic_pattern("C3H6O";R=1e9,print=false)
        @test length(actual) == length(expected)
        for (observed, ref) in zip(actual, expected)
            @test observed[1] ≈ ref[1] atol=1e-10 rtol=0
            @test observed[2] ≈ ref[2] atol=1e-11 rtol=0
        end
    end

    @testset "Principal isotope masses" begin
        for (formula, expected) in (("Fe",55.93493633), ("Li",7.0160034366), ("B",11.00930536),
                                    ("[54Fe]",53.93960899), ("[13C]C2H6O",59.04521964802),
                                    ("C3H6O",58.04186481295), ("C1000",12000.0))
            @test monoisotopic_mass(formula) ≈ expected atol=1e-10 rtol=0
        end
        @test monoisotopic_mass_protonated("Fe") ≈ 55.93493633+1.00782503223-0.0005485 atol=1e-12
        @test monoisotopic_mass("C3D6O") == monoisotopic_mass("C3[2H]6O")
        for invalid in ("", "2H", "[999C]", "Xx", "H0", "(H")
            @test_throws ArgumentError monoisotopic_mass(invalid)
        end
        for (element, mass) in (("Fe",55.93493633), ("Li",7.0160034366), ("B",11.00930536))
            matches = find_formula(mass;atom_pool=Dict(element=>1),tolerance_ppm=1e-6)
            @test any(c -> c.formula == element && isapprox(c.mz,mass;atol=1e-12), matches)
        end
        adduct_mass = 58.04186481295 + 7.0160034366 - 0.0005485
        matches = find_formula(adduct_mass;ion="[M+Li]+",atom_pool=Dict("C"=>3,"H"=>6,"O"=>1),tolerance_ppm=1e-6)
        @test any(c -> c.formula == "C3H6O" && isapprox(c.mz,adduct_mass;atol=1e-12), matches)
    end

    @testset "Resolution and centroid conservation" begin
        group = IsotopicCalc.IsotopicPatternModule.group_by_resolution
        peaks = [(100.0,1.0),(100.005,0.5),(100.02,0.25)]
        merged = group(reverse(peaks),10000)
        @test length(merged) == 2
        @test merged[1][1] ≈ (100.0+0.5*100.005)/1.5 atol=1e-12
        @test merged[1][2] == 1.5
        @test sum(last, merged) == sum(last, peaks)
        @test sum(m*a for (m,a) in merged) ≈ sum(m*a for (m,a) in peaks)
        @test group(peaks,Inf) == peaks
        @test isempty(group(Tuple{Float64,Float64}[],10000))
        @test length(group([(100.0,1.0),(100.009,1.0),(100.018,1.0)],10000)) == 2

        pattern = isotopic_pattern("C2H4";R=1000,print=false)
        m1 = filter(p -> 29 < p[1] < 30, pattern)
        @test length(m1) == 1
        carbon = 2*0.0107/0.9893
        hydrogen = 4*0.000115/0.999885
        @test m1[1][2] ≈ carbon+hydrogen atol=1e-12
        base = 24+4*1.00782503223
        @test m1[1][1] ≈ base+(carbon*1.00335483507+hydrogen*1.00627674589)/(carbon+hydrogen) atol=1e-10
        @test isotopic_pattern("C2H4";R=1,abundance_cutoff=1,print=false)[1][2] == 1.0
    end

    @testset "Adduct loss availability" begin
        invalid_mass = 12+2*15.99491461957-1.00782503223+0.0005485
        @test isempty(find_formula(invalid_mass;ion="[M-H]-",atom_pool=Dict("C"=>1,"O"=>2),tolerance_ppm=1))
        @test_throws ArgumentError isotopic_pattern("CO2";ion="[M-H]-",print=false)
        water_loss = 1.00782503223+15.99491461957+0.0005485
        matches = find_formula(water_loss;ion="[M-H]-",atom_pool=Dict("H"=>2,"O"=>1),tolerance_ppm=1e-6)
        @test any(c -> c.formula == "H2O" && c.charge == -1, matches)
        @test isotopic_pattern("H2O";ion="[M-H]-",print=false)[1][1] ≈ water_loss atol=1e-12
        @test isempty(find_formula(15.99491461957-1.00782503223+0.0005485;
                     ion="[M-3H]-",atom_pool=Dict("H"=>2,"O"=>1),tolerance_ppm=1))
        @test isempty(find_formula(1.00782503223-0.0005485;ion="[M+H]+",atom_pool=Dict("C"=>0),tolerance_ppm=1))
    end

    @testset "Valence and optional empirical limits" begin
        cf4_mass = 12+4*18.99840316273
        @test any(c -> c.formula == "CF4",find_formula(cf4_mass;atom_pool=Dict("C"=>1,"F"=>4),tolerance_ppm=1e-6))
        impossible = 24+6*1.00782503223+2*18.99840316273
        for heuristics in (true,false)
            @test isempty(find_formula(impossible;atom_pool=Dict("C"=>2,"H"=>6,"F"=>2),tolerance_ppm=1,heuristics=heuristics))
        end
        # Carbon disulfide is valid but lies outside the default S/C heuristic.
        cs2 = 12+2*31.9720711744
        @test isempty(find_formula(cs2;atom_pool=Dict("C"=>1,"S"=>2),tolerance_ppm=1))
        @test any(c -> c.formula == "CS2",find_formula(cs2;atom_pool=Dict("C"=>1,"S"=>2),tolerance_ppm=1,heuristics=false))
        # Isotope labels must not change neutral valence/parity.
        @test any(c -> c.formula == "CH3D", find_formula(12+3*1.00782503223+2.01410177812;
                  atom_pool=Dict("C"=>1,"H"=>3,"D"=>1),tolerance_ppm=1))
    end

    @testset "Cutoff boundaries and input validation" begin
        for printing in (false,true)
            zero = isotopic_pattern("C";abundance_cutoff=0,print=printing)
            @test length(zero) == 2
            @test zero[2][2] ≈ 0.0107/0.9893 atol=1e-14
            @test isotopic_pattern("C";abundance_cutoff=1,print=printing) == [(12.0,1.0)]
        end
        # A rare atomic isotope can contribute to the most abundant molecule.
        @test length(isotopic_pattern("Br2";abundance_cutoff=1,print=false)) == 1
        for cutoff in (-1,1.1,Inf,NaN)
            @test_throws ArgumentError isotopic_pattern("C";abundance_cutoff=cutoff,print=false)
        end
        for resolution in (-1,0,NaN)
            @test_throws ArgumentError isotopic_pattern("C";R=resolution,print=false)
        end
        for mz in (0,-1,NaN,Inf)
            @test_throws ArgumentError find_formula(mz)
        end
        for tolerance in (-1,NaN,Inf)
            @test_throws ArgumentError find_formula(12;tolerance_ppm=tolerance)
        end
        @test_throws ArgumentError find_formula(12;atom_pool=Dict("C"=>-1))
        @test_throws ArgumentError find_formula(12;atom_pool=Dict("Xx"=>1))
        @test any(c -> c.formula == "C2",find_formula(23.9999;atom_pool=Dict("C"=>2),tolerance_ppm=5))
    end
end
