using ACEbonds
using Test

@testset "ACEbonds.jl" begin
    # Write your tests here.
    @testset "Bonds basics" begin include("test_bonds.jl"); end 
    @testset "Invariant Cylindrical" begin include("test_invcyl.jl"); end
end
