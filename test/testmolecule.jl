function test_molecule()
    mol = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond1 = Bond(2)
    bond2 = Bond(1)
    add_bond!(mol, 1, 2, bond1)
    add_bond!(mol, 1, 3, bond2)
    return mol
end

function test_natom01()
    mol = test_molecule()
    n_atom = natom(mol)
    @test n_atom == 3
end

function test_nbond01()
    mol = test_molecule()
    n_bond = nbond(mol)
    @test n_bond == 2
end

function test_moltosmiles01()
    ssmiles = [
        "CCC",
        "C=CC",
        "CC=C",
        "C1CCCCC1",
        "C1CC(=CCC)CCC1"
    ]
    for smiles₀ ∈ ssmiles
        mol = smilestomol(smiles₀)
        smiles₁ = moltosmiles(mol)
        @test smiles₀ == smiles₁
    end
end

@testset "TransMol" begin
    test_natom01()
    test_nbond01()
    test_moltosmiles01()
end