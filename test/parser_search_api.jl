@testset "Token parser and integration API" begin
    P = IsotopicCalc.IsotopicPatternModule
    @testset "Token boundaries and round trips" begin
        for bad in ("C[13C]l", "C(H)l", "C[13C]a", "[13C]h", "C[]H", "()", "(())", ")C(",
                    "C(H", "C)H", "C 2", "C²", "C[13C]2l", "C()H", "[C13]", "[13(C)]", "[13C](H)0")
            @test_throws ArgumentError monoisotopic_mass(bad)
        end
        for formula in ("[13C]C2H6O", "([13C]H2)3", "C3D6O", "([2H]2O)4", "[37Cl]Cl", "Ca(OH)2", "((CH3)2)2")
            parsed = P.parse_formula(formula)
            formatted = P.sum_formula(parsed)
            @test P.parse_formula(formatted) == parsed
            @test monoisotopic_mass(formatted) == monoisotopic_mass(formula)
        end
        @test P.sum_formula(P.parse_formula("[13C]C2H6O")) == "C2H6O[13C]"
        @test P.parse_formula("(C[13C]H2)3O") == Dict("C"=>3,"13C"=>3,"H"=>6,"O"=>1)
        @test P.parse_formula("(CH2)1000000000") == Dict("C"=>1000000000,"H"=>2000000000)
        @test P.parse_formula(repeat("(",1000)*"H"*repeat(")",1000)) == Dict("H"=>1)
        maximum_count = string(typemax(Int))
        @test_throws ArgumentError P.parse_formula("(H$maximum_count)2")
        @test_throws ArgumentError P.parse_formula("H$(maximum_count)H")
        @test_throws ArgumentError P.parse_formula("H$(maximum_count)0")
    end

    @testset "Strings, output destinations and protonated keywords" begin
        formula = SubString("!C3H6O!",2,6)
        ion = SubString("![M+H]+!",2,7)
        @test monoisotopic_mass(formula) == monoisotopic_mass("C3H6O")
        @test monoisotopic_mass_protonated(formula) == monoisotopic_mass_protonated("C3H6O")
        @test isotopic_pattern(formula;ion=ion,print=false) == isotopic_pattern("C3H6O";ion="[M+H]+",print=false)
        output = IOBuffer()
        result = isotopic_pattern_protonated(formula;R=Inf,abundance_cutoff=0.01,io=output)
        @test result == isotopic_pattern("C3H6O";ion="[M+H]+",R=Inf,abundance_cutoff=0.01,print=false)
        @test occursin("C3H6O.[M+H]+",String(take!(output)))
        @test all(p -> p[2] >= 0.01,result)
        isotopic_pattern_protonated(formula;io=output,print=false)
        @test isempty(take!(output))
        @test_throws ArgumentError isotopic_pattern_protonated(formula;ion="[M+Na]+")
        isotopic_pattern("[13C]C2H6O";io=output)
        @test occursin("C2H6O[13C]",String(take!(output)))
        matches = find_formula(59.0491;ion=ion,io=output)
        @test occursin("Matching formulas:",String(take!(output)))
        quiet = find_formula(59.0491;ion=ion,print=false,io=output)
        @test [(c.formula,c.mz,c.ppm) for c in quiet] == [(c.formula,c.mz,c.ppm) for c in matches]
        @test isempty(take!(output))
        # A separate stdout capture verifies that an explicit IO destination is
        # respected, rather than merely checking the supplied buffer is empty.
        mktemp() do path, stream
            redirect_stdout(stream) do
                find_formula(59.0491;ion=ion,print=false)
                find_formula(59.0491;ion=ion,io=output)
                isotopic_pattern_protonated(formula;print=false)
                isotopic_pattern(formula;io=output)
            end
            flush(stream)
            @test isempty(read(path,String))
        end
    end
end

@testset "Mass-pruned search against exhaustive enumeration" begin
    # Independent exhaustive reference: no branch bounds, string construction,
    # isotope lookups, or production search helpers are used to decide matches.
    pool = Dict("C"=>4,"H"=>10,"Cl"=>2,"N"=>2,"O"=>2)
    records = Tuple{NTuple{5,Int},Float64}[]
    for c in 0:4, h in 0:10, cl in 0:2, n in 0:2, o in 0:2
        counts = (c,h,cl,n,o)
        all(iszero,counts) && continue
        dbe2 = 2c+2+n-h-cl
        (dbe2 < 0 || isodd(dbe2)) && continue
        # This is the canonical element order used for final exact masses.
        mass = sum((12.0c,34.968852682cl,1.00782503223h,14.00307400443n,15.99491461957o))
        push!(records,(counts,mass))
    end
    scenarios = (("[M]",0,0.0,0,0), ("[M+H]+",1,1.00782503223,0,0),
                 ("[M+2H]2+",2,2*1.00782503223,0,0), ("[M]3+",3,0.0,0,0),
                 ("[M-H]-",-1,-1.00782503223,1,0), ("[M-Cl]-",-1,-34.968852682,0,1))
    # Several target masses with matches, both edges, and points outside 5 ppm.
    for (ion,charge,delta,min_h,min_cl) in scenarios
        available = filter(r -> r[1][2] >= min_h && r[1][3] >= min_cl, records)
        target_indices = unique(round.(Int,range(1,length(available);length=8)))
        for index in target_indices, offset_ppm in (-6.0,-5.0,0.0,5.0,6.0)
            theoretical = (available[index][2]+delta-charge*0.0005485)/max(abs(charge),1)
            observed = theoretical*(1+offset_ppm/1e6)
            observed > 0 || continue
            expected = Set{NTuple{5,Int}}()
            for (counts,mass) in available
                mz = (mass+delta-charge*0.0005485)/max(abs(charge),1)
                mz > 0 && abs((observed-mz)/mz*1e6) <= 5 && push!(expected,counts)
            end
            matches = find_formula(observed;ion=ion,atom_pool=pool,tolerance_ppm=5,heuristics=false,print=false)
            actual = Set(Tuple(get(IsotopicCalc.IsotopicPatternModule.parse_formula(item.formula),e,0)
                               for e in ("C","H","Cl","N","O")) for item in matches)
            @test actual == expected
            @test issorted(abs.([item.ppm for item in matches]))
        end
    end
    for tolerance in (0.0, 1e6, 2e6)
        target = 24.0
        expected = Set(counts for (counts,mass) in records if abs((target-mass)/mass*1e6) <= tolerance)
        matches = find_formula(target;atom_pool=pool,tolerance_ppm=tolerance,heuristics=false,print=false)
        actual = Set(Tuple(get(IsotopicCalc.IsotopicPatternModule.parse_formula(item.formula),e,0)
                           for e in ("C","H","Cl","N","O")) for item in matches)
        @test actual == expected
    end
    @test isempty(find_formula(12;atom_pool=Dict{String,Int}(),print=false))
    @test isempty(find_formula(12;atom_pool=Dict("H"=>0),print=false))
    @test_throws ArgumentError find_formula(1;ion="[M+Na]+",atom_pool=Dict("Xx"=>1),print=false)
end
