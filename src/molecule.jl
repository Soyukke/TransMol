export Molecule
export add_atom!, add_atom, add_bond, add_bond!
export natom, nbond, writesdf, smilestomol
using Printf
using LightGraphs, MetaGraphs 
using MolecularGraph: sdftomol, sdfilewriter, coordgen, drawsvg

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

function get_atom(m::Molecule, i::Integer)
    atom = get_prop(m.graph, i, :atom)
    return atom
end

function get_coords(m::Molecule, i::Integer)
    if !has_prop(m.graph, i, :coords)
        return Float64[0, 0, 0]
    else
        return get_prop(m.graph, i, :coords)
    end
end

function set_coords(m::Molecule, i::Integer, coords::Vector{T}) where T <: Real
    set_prop!(m.graph, i, :coords, coords)
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
全原子
"""
function atoms(mol::Molecule)
    natom_ = natom(mol)
    allatoms = map(i -> get_prop(mol.graph, i, :atom), 1:natom_)
    # 原子番号でsort
    sort!(allatoms, by = atom -> atom.number)
    return allatoms
end

"""
原子数
"""
function natom(mol::Molecule)
    return nv(mol.graph)
end

"""
結合数
"""
function nbond(mol::Molecule)
    return ne(mol.graph)
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
    mols = Molecule[]
    n_atom = natom(mol)
    for i ∈ 1:n_atom, j ∈ 1:n_atom
        nf1 = nfree(mol, i)
        nf2 = nfree(mol, j)
        order = bondorder(mol, i, j)
        nfreemin = min(nf1, nf2)
        bonds = filter(x->x ≤ nfreemin, bonddict[order])
        for transorder ∈ bonds
            newbond = Bond(transorder)
            push!(mols, add_bond(mol, i, j, newbond))
        end
        # @show order, nfreemin, filter(x->x ≤ nfreemin, bonddict[order])
    end
    return mols
end

"""
結合削除
singlebond => [None]                        : 結合次数1
doublebond => [None, singlebond]            : 結合次数2
triplebond => [None, singlebond, doublebond]: 結合次数3
"""
function removable_bond(mol::Molecule)
    # 結合削除dict
    bonddict = Dict(0 => [], 1 => [0], 2 => [0, 1], 3 => [0, 1, 2])
    n_atom = natom(mol)
    mols = Molecule[]
    list = []
    for i ∈ 1:n_atom, j ∈ 1:n_atom
        order = bondorder(mol, i, j)
        push!(list, bonddict[order])
    end
    return list
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
    @show removable_bond(mol)
    mol
end

function example2()
    mol = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond1 = Bond(2)
    bond2 = Bond(1)
    add_bond!(mol, 1, 2, bond1)
    add_bond!(mol, 1, 3, bond2)
    fpath = joinpath(@__DIR__, "..", "devtest", "mol.sdf")
    writesdf(fpath, mol)
    mol2 = sdftomol(fpath)
    svglines = drawsvg(mol2, 600, 600)
    fsvgpath = joinpath(@__DIR__, "..", "devtest", "mol.svg")
    write(fsvgpath, svglines)
end

"""
Molecule -> sdf
https://en.wikipedia.org/wiki/Chemical_table_file
"""
function sdf(mol::Molecule)
    name = "mol name"
    programname = "program name"
    comment = "hoge"
    n_atom = natom(mol)
    n_bond = nbond(mol)
    natomlist = 0
    chiralfrag = 0
    nstextentries = 0
    # M ENDはカウントされる
    nprops = 1
    molversion = "V2000"
    numbers = [n_atom, n_bond, natomlist, chiralfrag, nstextentries, nprops]
    extra = " 0  0  0  0999"
    strs = map(x -> @sprintf(" %d ", x), numbers)
    counts = " " * join(strs, "") * extra * " $(molversion)"
    atom_block = String[]
    for i ∈ 1:n_atom
        atom = get_atom(mol, i)
        x, y, z = get_coords(mol, i)
        name = "$(string(atom.name)) "
        massdiff = 0
        charge = 0
        sss = 0
        nhydrogen = 0
        bbb = 0
        valence = 0
        numbers = [massdiff, charge, sss, nhydrogen, bbb, valence, 0, 0, 0, 0, 0, 0,]
        nstrs = map(x -> @sprintf(" %d", x), numbers)
        line = "" *
        @sprintf("% 10.4f", x) *
        @sprintf("% 10.4f", y) *
        @sprintf("% 10.4f", z) *
        " " *
        name * " " * join(nstrs, " ")
        push!(atom_block, line)
        # hhh, rrr, iii, mmm, nnn, eee = 0, 0, 0, 0, 0, 0
    end

    # 結合テーブル
    bond_block = String[]
    for i ∈ 1:n_atom, j ∈ 1:n_atom
        nf1 = nfree(mol, i)
        nf2 = nfree(mol, j)
        order = bondorder(mol, i, j)
        if order == 0
            continue
        end
        numbers = [i, j, order, 0, 0, 0, 0]

        strs = map(x -> @sprintf(" %d ", x), numbers)
        line = " " * join(strs, "")
        push!(bond_block, line)
        # @show order, nfreemin, filter(x->x ≤ nfreemin, bonddict[order])
    end
    prop_block = String[]
    mend = "M  END"
    push!(prop_block, mend)
    lines = [
        name, programname, comment, counts,
        atom_block..., bond_block..., prop_block...
    ]
    return lines
end

"""
座標をセットしてsdfファイルを出力する
"""
function sdfwithcoords(mol::Molecule)
    lines = sdf(mol)
    sdfmol = sdftomol(lines)
    # 描画用座標をセットする
    coords, styles = coordgen(sdfmol)
    for i ∈ 1:natom(mol)
        set_coords(mol, i, Float64[coords[i, :]..., 0])
    end
    lines = sdf(mol)
    text = join(lines, "\n")
end

function writesdf(filename::AbstractString, mol::Molecule)
    text = sdfwithcoords(mol)
    write(filename, text)
end

"""
smiles(mol::Molecule)::String

SMILES文字列
"""
function smiles(mol::Molecule)
end

atomdict = Dict()
for a in atomdata
    push!(atomdict, string(a[begin]) => a)
end
bonddict = Dict(""=>1, "="=>2, "#"=>3)

function smilestomol(smiles::String)
    mol = Molecule()
    x = split(smiles, "")
    next = iterate(x)
    smilestomol(mol, x, 1)
end

"""
再帰的にMoleculeを作る
"""
function smilestomol(mol::Molecule, x::Vector, i)
    next = iterate(x, i)
    while next !== nothing 
        t, i = next
        # nestを抜ける
        if t == ")"
            println("out next")
            return mol, t
        elseif t == "("
            # 結合をつなげる原子indexを保持する？
            println(t, "in nest")
            smilestomol(mol, x, i)
        else t ∈ ["1", "2", "3", "4", "5", "6", "7", "8"]
        end
        # 原子と一致
        if haskey(atomdict, t)
            println("原子")
            atomidx += 1
        elseif haskey(bonddict, t)
            println("結合")
        end
        next = iterate(x, i)
    end
    return mol, nothing
end