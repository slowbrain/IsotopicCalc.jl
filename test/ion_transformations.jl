@testset "Hydrogen subtraction independent of charge" begin
    # Independent masses from the bundled principal isotopes, without ion helpers.
    h, o, electron = 1.00782503223, 15.99491461957, 0.0005485
    precursor_mass = 36 + 6h + o
    pool = Dict("C"=>3, "H"=>6, "O"=>1)
    for (notation, removed, charge) in (("[M-H]-", 1, -1), ("[M-H]+", 1, 1),
                                       ("[M-H]", 1, 0), ("[M-2H]2-", 2, -2))
        expected = (precursor_mass - removed*h - charge*electron) / max(abs(charge), 1)
        spec = IonSpec(remove=removed == 1 ? "H" : "H2", charge=charge)
        @test ion_mass("C3H6O"; ion=notation) ≈ expected atol=1e-12
        @test ion_mass("C3H6O"; ion=spec) ≈ expected atol=1e-12
        pattern = isotopic_pattern("C3H6O"; ion=notation, R=Inf, annotated=true, print=false)
        @test first(pattern.peaks).mass ≈ expected atol=1e-12
        @test pattern.ion_composition == Dict("C"=>3, "H"=>6-removed, "O"=>1)
        @test pattern.charge == charge
        @test pattern.settings.ion.losses == Dict("H"=>removed)

        search = find_formula(expected; ion=notation, atom_pool=pool,
                              tolerance_ppm=1e-6, diagnostics=true, print=false)
        compound = only(search.matches)
        @test compound.formula == "C3H6O"
        @test compound.mz ≈ expected atol=1e-12
        @test compound.charge == charge
        @test ion_mass(compound.formula; ion=compound.ion) ≈ expected atol=1e-12
        @test search.diagnostics.status == :matched

        # Ranking must apply the same loss to isotope predictions as to search.
        observed = [(peak.mass, peak.abundance) for peak in pattern.peaks]
        ranked = only(rank_formulas(expected, observed;
            search_options=(ion=notation, atom_pool=pool, tolerance_ppm=1e-6),
            pattern_options=(R=Inf,)))
        @test ranked.compound.formula == "C3H6O"
        @test ranked.isotope_similarity ≈ 1 atol=1e-12
        @test isempty(ranked.missing_peaks)
        @test isempty(ranked.unexplained_peaks)

        exported = result_dict(ranked)
        @test exported["schema_version"] == 2
        @test exported["data"]["compound"]["ion"] == compound.ion
        @test !haskey(exported["data"]["compound"], "adduct")
        @test exported["data"]["pattern"]["settings"]["ion"]["losses"] == Dict("H"=>removed)

        @test_throws ArgumentError ion_mass("CO2"; ion=notation)
        @test_throws ArgumentError isotopic_pattern("CO2"; ion=notation, print=false)
        infeasible = find_formula(40; ion=notation, atom_pool=Dict("C"=>4),
                                  diagnostics=true, print=false)
        @test isempty(infeasible.matches)
        @test infeasible.diagnostics.status == :infeasible_ion
    end

    # Labeled hydrogen is not available for an unlabeled hydrogen subtraction.
    @test_throws ArgumentError ion_mass("C3[2H]6O"; ion="[M-H]+")
    @test_throws ArgumentError isotopic_pattern("C3[2H]6O"; ion="[M-H]-", print=false)
    @test ion_mass("C3[2H]6O"; ion="[M-[2H]]+") ≈ 36 + 5*2.01410177812 + o - electron atol=1e-12
    @test_throws ArgumentError ion_mass("H"; ion="[M-2H]-")
    @test_throws ArgumentError ion_mass("H"; ion="[M-H]") # Empty product has no positive mass.

    # Both forms use the same full label in user-facing output.
    output = IOBuffer()
    isotopic_pattern("C3H6O"; ion=IonSpec(remove="H", charge=1), io=output)
    @test occursin("C3H6O.[M-H]+", String(take!(output)))
    find_formula(precursor_mass-h-electron; ion="[M-H]+", atom_pool=pool, io=output)
    @test occursin("[M-H]+", String(take!(output)))
end

@testset "Explicit notation and removal of legacy API" begin
    @test ion_mass("H2O") == ion_mass("H2O"; ion="[M]")
    @test ion_mass("H2O"; ion=IonSpec()) == monoisotopic_mass("H2O")
    @test string(IonSpec()) == "[M]"
    @test string(IonSpec(remove="H", charge=1)) == "[M-H]+"
    @test string(IonSpec(remove="H", charge=-1)) == "[M-H]-"
    @test string(IonSpec(remove="H")) == "[M-H]"
    @test string(IonSpec(add="Cl", charge=-1)) == "[M+Cl]-"
    for notation in ("[M]", "[M-H]+", "[M-2H]2-", "[M+NH4]+", "[M+H-H2O]+",
                     "[M+[13C]-[2H]]+", "[M]2+", "[M]2-", "[M+H-H2O]")
        spec = IonSpec(notation)
        parsed = IonSpec(string(spec))
        @test parsed.additions == spec.additions
        @test parsed.losses == spec.losses
        @test parsed.charge == spec.charge
    end
    for legacy in ("", "H+", "H-", "Cl-", "NH4+", "2H+2", "2H-2", "+", "-", "+2", "-2", "M-H")
        @test_throws ArgumentError IonSpec(legacy)
        @test_throws ArgumentError ion_mass("H2O"; ion=legacy)
        @test_throws ArgumentError isotopic_pattern("H2O"; ion=legacy, print=false)
        @test_throws ArgumentError find_formula(18; ion=legacy, print=false)
    end
    @test_throws MethodError ion_mass("H2O"; adduct="[M-H]-")
    @test_throws MethodError isotopic_pattern("H2O"; adduct="[M-H]-", print=false)
    @test_throws MethodError find_formula(18; adduct="[M-H]-", print=false)
    @test_throws MethodError isotopic_pattern_protonated("H2O"; adduct="H+")
    for name in (:AdductInfo, :parse_adduct, :get_adduct_info)
        @test !isdefined(IsotopicCalc.FindFormulaModule, name)
    end
    for name in (:PatternAdductInfo, :parse_pattern_adduct)
        @test !isdefined(IsotopicCalc.IsotopicPatternModule, name)
    end
end
