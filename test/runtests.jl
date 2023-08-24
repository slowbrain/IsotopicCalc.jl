using IsotopicCalc
using Test

@testset "IsotopicCalc.jl" begin
    # Test the formula search
    mz_input = 33.034
    tolerance = 20.0  # Tolerance in ppm
    atom_pool = Dict("C" => 5, "H" => 10, "N" => 2, "O" => 3)
    adduct = "M+H"
    charge = 1

    @test find_compounds(mz_input, tolerance, atom_pool, adduct, charge) == [Compound("CH4O", "M+H", 1, 33.03349128072001, 15.400106385027454)]
end
