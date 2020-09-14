export Molecule
export add_atom!
using LightGraphs, MetaGraphs

struct Molecule
    graph::MetaGraph
end

function Molecule()
    g = MetaGraph()
    return Molecule(g)
end

function Base.copy(mol::Molecule)
    return Molecule(copy(mol.graph))
end

function Base.deepcopy(mol::Molecule)
    return Molecule(deepcopy(mol.graph))
end

function Base.isequal(m1::Molecule, m2::Molecule)
    return m1.graph == m2.graph
end

"""
atomをmoleculeに追加する
Bool値をreturnする
`true`の場合追加頂点追加成功
"""
function add_atom!(m::Molecule, a::Atom)
    dict = Dict(:atom=>a)
    add_vertex!(m.graph, dict)
end

function add_atom(m::Molecule, a::Atom)
    m2 = deepcopy(m)
    add_atom!(m2, a)
    return m2
end

"""
分子`m`へ原子`i`と原子`j`に結合`b`を追加する
"""
function add_bond!(m::Molecule, i::Integer, j::Integer, b::Bond)
    dict = Dict(:bond=>b)
    add_edge!(m.graph, i, j, dict)
end

function add_bond(m::Molecule, i::Integer, j::Integer, b::Bond)
    m2 = deepcopy(m)
    add_bond!(m2, i, j, b)
    return m2
end


"""
atom `i` に関する残り原子価数
atom `i`に接続されているbondsのbond order合計をatom `i`の原子価数から差し引いた値
"""
function nfree(m::Molecule, i::Integer)
    # iに隣接する原子すべてとの結合字数の総和を計算する
    indices = neighbors(m.graph, i)
    n = reduce(indices) do x, j
        bond = get_prop(m.graph, i, j, :bond)
        y = bond.order
        return x + y
    end
    atom = get_prop(m.graph, i, :atom)
    valence = atom.valence
    return valence - n
end

function bondorder(m::Molecule, i::Integer, j::Integer)
    try
        bond = get_prop(m.graph, i, j, :bond)
        if bond !== nothing
            return bond.order
        end
    catch e
        return Int(0)
    end
    # finallyでreturn すると必ずfinallyでreturnされるのでこうする
end

"""
原子数
"""
function natom(mol::Molecule)
    return nv(mol.graph)
end

"""
原子追加
すべての原子a ∈ molについて水素を除いた残価数を計算
残価数 `nfree` が1以上なら，原子b ∈ 原子リストが追加可能である．
また，原子追加後に同じ分子は削除してuniqueな原子追加を残す
"""
function addable_atom(mol::Molecule)
    list = []
    atoms = vertices(mol.graph)
    for (i, a) ∈ enumerate(atoms)
        nf = nfree(mol, i)
        if 1 ≤ nf
            d = [i, atomlist]
            push!(list, d)
        end
    end
    println(list)
end

"""
結合追加
None => [singlebond, doublebond, triplebond] : 原子間結合なし
singlebond => [doublebond, triplebond]       : 結合次数1
doublebond => [triplebond]                   : 結合次数2
原子i, 原子jについてそれぞれ残原子価が1以上のときに結合追加可能なパターンを列挙する
結合追加後の分子が同じものは取り除く
"""
function addable_bond(mol::Molecule)
    bonddict = Dict(0 => [1, 2, 3], 1 => [2, 3], 2 => [3])
    n_atom = natom(mol)
    for i ∈ 1:n_atom, j ∈ 1:n_atom
        nf1 = nfree(mol, i)
        nf2 = nfree(mol, j)
        order = bondorder(mol, i, j)
        nfreemin = min(nf1, nf2)
        @show order, nfreemin, filter(x->x ≤ nfreemin, bonddict[order])
    end
end

"""
結合削除
"""
function removable_bond(mol::Molecule)
end

function example1()
    mol = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond1 = Bond(2)
    bond2 = Bond(1)
    add_bond!(mol, 1, 2, bond1)
    add_bond!(mol, 1, 3, bond2)
    nf = nfree(mol, 1)
    @show nf
    @show natom(mol)
    @show bondorder(mol, 1, 2)
    @show addable_atom(mol)
    @show addable_bond(mol)
    mol
end