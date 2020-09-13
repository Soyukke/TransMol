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

"""
atomをmoleculeに追加する
Bool値をreturnする
`true`の場合追加頂点追加成功
"""
function add_atom!(m::Molecule, a::Atom)
    dict = Dict(:atom=>a)
    add_vertex!(m.graph, dict)
end

"""
分子`m`へ原子`i`と原子`j`に結合`b`を追加する
"""
function add_bond!(m::Molecule, i::Integer, j::Integer, b::Bond)
    dict = Dict(:bond=>b)
    add_edge!(m.graph, i, j, dict)
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

function example1()
    mol = Molecule()
    C = Atom(6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond = Bond(2)
    add_bond!(mol, 1, 2, bond)
    nf = nfree(mol, 1)
    @show nf
    mol
end