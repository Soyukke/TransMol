import TransMol:aromatic2atom

function test_aromatic2atom()
    @test "Br" == aromatic2atom("br")
    @test "C" == aromatic2atom("c")
end

@testset "smiles.jl" begin
    test_aromatic2atom()
end