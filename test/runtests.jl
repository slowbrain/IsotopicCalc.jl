using IsotopicCalc
using Test

@testset "IsotopicCalc.jl" begin

    @testset "monoisotopicMass - Basic Formulas" begin
        # Test simple molecules (tolerance accounts for R=10000 resolution rounding)
        @test isapprox(monoisotopicMass("H2O"), 18.010565, atol=1e-4)
        @test isapprox(monoisotopicMass("CO2"), 43.989829, atol=1e-4)
        @test isapprox(monoisotopicMass("CH4"), 16.031300, atol=1e-4)
        @test isapprox(monoisotopicMass("C3H6O"), 58.041865, atol=1e-4)  # Acetone
        @test isapprox(monoisotopicMass("C6H12O6"), 180.063388, atol=1e-3)  # Glucose

        # Test single atoms
        @test isapprox(monoisotopicMass("C"), 12.0, atol=1e-4)
        @test isapprox(monoisotopicMass("H"), 1.007825, atol=1e-5)
        @test isapprox(monoisotopicMass("O"), 15.994915, atol=1e-4)
        @test isapprox(monoisotopicMass("N"), 14.003074, atol=1e-4)
    end

    @testset "monoisotopicMass - Parentheses Expansion" begin
        # Test equivalent formulas with parentheses
        @test isapprox(monoisotopicMass("(CH3)2CO"), monoisotopicMass("C3H6O"), atol=1e-6)
        @test isapprox(monoisotopicMass("Ca(OH)2"), monoisotopicMass("CaO2H2"), atol=1e-6)
        @test isapprox(monoisotopicMass("CH3COCH3"), monoisotopicMass("C3H6O"), atol=1e-6)

        # Test nested parentheses
        @test isapprox(monoisotopicMass("((CH3)2)2"), monoisotopicMass("C4H12"), atol=1e-6)

        # Test multiple groups
        @test isapprox(monoisotopicMass("(CH2)3"), monoisotopicMass("C3H6"), atol=1e-6)
    end

    @testset "monoisotopicMass - Isotope Specifications" begin
        # Test specific isotopes
        c12_mass = monoisotopicMass("C")
        c13_mass = monoisotopicMass("[13C]")
        @test c13_mass > c12_mass  # 13C should be heavier
        @test isapprox(c13_mass - c12_mass, 1.003355, atol=1e-4)  # Mass difference

        # Test deuterium
        h_mass = monoisotopicMass("H")
        d_mass = monoisotopicMass("[2H]")
        @test d_mass > h_mass
        @test isapprox(d_mass - h_mass, 1.006277, atol=1e-4)

        # Test D shorthand for deuterium
        @test isapprox(monoisotopicMass("D"), monoisotopicMass("[2H]"), atol=1e-6)
        @test isapprox(monoisotopicMass("C3D6O"), monoisotopicMass("C3[2H]6O"), atol=1e-6)

        # Test mixed isotopes
        acetone_normal = monoisotopicMass("C3H6O")
        acetone_13C = monoisotopicMass("[13C]C2H6O")
        @test isapprox(acetone_13C - acetone_normal, 1.003355, atol=1e-4)
    end

    @testset "monoisotopicMassProtonated" begin
        # Test protonated masses
        acetone = monoisotopicMass("C3H6O")
        acetone_protonated = monoisotopicMassProtonated("C3H6O")
        @test acetone_protonated > acetone
        @test isapprox(acetone_protonated - acetone, 1.007825, atol=1e-4)  # Proton mass

        # Compare with direct isotopicPattern call
        pattern_h = isotopicPattern("C3H6O"; adduct="H+", print=false)
        @test isapprox(monoisotopicMassProtonated("C3H6O"), pattern_h[1][1], atol=1e-6)
    end

    @testset "isotopicPattern - Basic Functionality" begin
        # Test that function returns expected data structure
        pattern = isotopicPattern("CH4"; print=false)
        @test isa(pattern, Vector{Tuple{Float64, Float64}})
        @test length(pattern) >= 1

        # First peak should be monoisotopic mass
        @test isapprox(pattern[1][1], monoisotopicMass("CH4"), atol=1e-6)

        # Most abundant peak should have abundance 1.0
        abundances = [a for (m, a) in pattern]
        @test maximum(abundances) ≈ 1.0

        # Masses should be sorted
        masses = [m for (m, a) in pattern]
        @test issorted(masses)
    end

    @testset "isotopicPattern - Adducts" begin
        # Test with different adducts
        base_mass = monoisotopicMass("C3H6O")

        # H+ adduct
        pattern_h = isotopicPattern("C3H6O"; adduct="H+", print=false)
        @test isapprox(pattern_h[1][1] - base_mass, 1.007825, atol=1e-4)

        # Na+ adduct
        pattern_na = isotopicPattern("C3H6O"; adduct="Na+", print=false)
        @test isapprox(pattern_na[1][1] - base_mass, 22.989769, atol=1e-4)

        # K+ adduct
        pattern_k = isotopicPattern("C3H6O"; adduct="K+", print=false)
        @test isapprox(pattern_k[1][1] - base_mass, 38.963706, atol=1e-4)

        # H- (deprotonated)
        pattern_minus = isotopicPattern("C3H6O"; adduct="H-", print=false)
        @test pattern_minus[1][1] < base_mass
    end

    @testset "isotopicPattern - Parameters" begin
        # Test abundance cutoff parameter
        pattern_default = isotopicPattern("C10H20O5"; print=false)
        pattern_strict = isotopicPattern("C10H20O5"; abundance_cutoff=1e-3, print=false)
        @test length(pattern_strict) <= length(pattern_default)

        # Test resolution parameter
        pattern_high_res = isotopicPattern("C3H6O"; R=100000, print=false)
        pattern_low_res = isotopicPattern("C3H6O"; R=1000, print=false)
        @test length(pattern_high_res) >= length(pattern_low_res)  # Higher R may resolve more peaks
    end

    @testset "isotopicPattern - Isotope Effects" begin
        # Molecules with more carbons should have more M+1 peak
        c3_pattern = isotopicPattern("C3H6O"; print=false)
        c10_pattern = isotopicPattern("C10H20O5"; print=false)

        # Find M+1 relative abundance (second peak typically)
        if length(c3_pattern) >= 2 && length(c10_pattern) >= 2
            c3_m1 = c3_pattern[2][2]
            c10_m1 = c10_pattern[2][2]
            @test c10_m1 > c3_m1  # More carbons = larger M+1 peak
        end
    end

    @testset "findFormula - Basic Functionality" begin
        # Test finding exact mass match
        mz = monoisotopicMass("C3H6O")
        matches = findFormula(mz; tolerance_ppm=10, atom_pool=Dict("C"=>5, "H"=>10, "O"=>3))

        @test length(matches) >= 1
        @test any(c -> c.formula == "C3H6O", matches)

        # Check that matches are sorted by ppm
        if length(matches) > 1
            ppms = [abs(c.ppm) for c in matches]
            @test issorted(ppms)
        end
    end

    @testset "findFormula - With Adducts" begin
        # Test the original formula search with adduct
        mz_input = 33.034
        tolerance = 20.0
        atom_pool = Dict("C" => 5, "H" => 10, "N" => 2, "O" => 3)
        adduct = "M+H"
        charge = 1

        results = findFormula(mz_input, tolerance_ppm=tolerance, atom_pool=atom_pool, adduct=adduct, charge=charge)
        @test length(results) >= 1
        @test any(c -> c.formula == "CH4O", results)

        # Check that found compound has correct fields
        ch4o = findfirst(c -> c.formula == "CH4O", results)
        if !isnothing(ch4o)
            compound = results[ch4o]
            @test compound.adduct == "M+H"
            @test compound.charge == 1
            @test abs(compound.ppm) <= tolerance
        end
    end

    @testset "findFormula - Tolerance Filtering" begin
        mz = 58.0419  # Acetone

        # Strict tolerance should return fewer results
        strict = findFormula(mz; tolerance_ppm=5, atom_pool=Dict("C"=>5, "H"=>10, "O"=>3, "N"=>3))
        loose = findFormula(mz; tolerance_ppm=1000, atom_pool=Dict("C"=>5, "H"=>10, "O"=>3, "N"=>3))

        @test length(strict) <= length(loose)
    end

    @testset "Compound Structure" begin
        # Test Compound struct creation and fields
        comp = Compound("C3H6O", "M+H", 1, 59.049141, -0.7)

        @test comp.formula == "C3H6O"
        @test comp.adduct == "M+H"
        @test comp.charge == 1
        @test comp.mz ≈ 59.049141
        @test comp.ppm ≈ -0.7

        # Test that Compound is properly exported
        @test isdefined(IsotopicCalc, :Compound)
    end

    @testset "Error Handling - Invalid Formulas" begin
        # Test invalid characters
        @test_throws ArgumentError isotopicPattern("C3H6@#")
        @test_throws ArgumentError isotopicPattern("ABC123!!")

        # Test unbalanced brackets
        @test_throws ArgumentError isotopicPattern("C3H6(O")
        @test_throws ArgumentError isotopicPattern("C3H6)O")
        @test_throws ArgumentError isotopicPattern("((CH3)")

        # Test invalid parameters
        @test_throws ArgumentError isotopicPattern("C3H6O"; abundance_cutoff=-0.1)
        @test_throws ArgumentError isotopicPattern("C3H6O"; R=-100)
        @test_throws ArgumentError isotopicPattern("C3H6O"; R=0)
    end

    @testset "Error Handling - Unknown Elements" begin
        # Test unknown element
        @test_throws ArgumentError isotopicPattern("Xx")  # Xx is not a real element
        @test_throws ArgumentError isotopicPattern("C3H6Zz")  # Zz is not real
    end

    @testset "Edge Cases - Single Atoms" begin
        # Test all common elements work as single atoms
        for element in ["H", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"]
            pattern = isotopicPattern(element; print=false)
            @test length(pattern) >= 1
            @test pattern[1][2] ≈ 1.0  # Most abundant peak
        end
    end

    @testset "Edge Cases - Large Molecules" begin
        # Test a moderately large molecule doesn't crash
        pattern = isotopicPattern("C20H30O10"; print=false)
        @test length(pattern) >= 1
        @test pattern[1][2] ≈ 1.0

        # Test protein-like formula
        pattern2 = isotopicPattern("C100H150N30O30S2"; abundance_cutoff=1e-3, print=false)
        @test length(pattern2) >= 1
    end

    @testset "Edge Cases - Elements with Multiple Digits" begin
        # Test formulas with 2-digit counts
        @test isapprox(monoisotopicMass("C10H22"), monoisotopicMass("C") * 10 + monoisotopicMass("H") * 22, atol=1e-4)
        @test isapprox(monoisotopicMass("C100"), monoisotopicMass("C") * 100, atol=1e-4)
    end

    @testset "Consistency Checks" begin
        # isotopicPattern and monoisotopicMass should agree
        formulas = ["H2O", "CO2", "C3H6O", "C6H12O6", "CH3COOH"]
        for formula in formulas
            pattern = isotopicPattern(formula; print=false)
            mass_direct = monoisotopicMass(formula)
            @test isapprox(pattern[1][1], mass_direct, atol=1e-6)
        end

        # Protonated functions should agree
        @test isapprox(
            isotopicPatternProtonated("C3H6O")[1][1],
            monoisotopicMassProtonated("C3H6O"),
            atol=1e-6
        )
    end

    @testset "Formula Equivalence" begin
        # Different representations of same molecule should give same mass
        acetone_forms = ["C3H6O", "(CH3)2CO", "CH3COCH3"]
        masses = [monoisotopicMass(f) for f in acetone_forms]
        @test all(m -> isapprox(m, masses[1], atol=1e-6), masses)
    end

end
