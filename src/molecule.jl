export Molecule
export add_atom!, add_atom, add_bond, add_bond!
export natom, nbond, writesdf, smilestomol
export atoms, bonds
export natom, nbond, writesdf
using Printf
using LightGraphs, MetaGraphs 
using MolecularGraph: sdftomol, sdfilewriter, coordgen, drawsvg, GraphMol, SmilesAtom, SmilesBond, graphmol
using Combinatorics

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

"""
i, jの結合を削除する
"""
function remove_bond!(m::Molecule, i::Integer, j::Integer)
    rem_edge!(m.graph, i, j)
end

function add_bond(m::Molecule, i::Integer, j::Integer, b::Bond)
    m2 = deepcopy(m)
    add_bond!(m2, i, j, b)
    return m2
end

function remove_bond(m::Molecule, i::Integer, j::Integer)
    m₂ = deepcopy(m)
    remove_bond!(m₂, i, j)
    return m₂
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
全結合
"""
function bonds(mol::Molecule)
    bondmaps = collect(edges(mol.graph))
    B = map(bondmaps) do bond
        get_prop(mol.graph, bond.src, bond.dst, :bond)
    end
    return B, bondmaps
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
    mols = Molecule[]
    atoms = vertices(mol.graph)
    for (i, a) ∈ enumerate(atoms)
        nf = nfree(mol, i)
        if 1 ≤ nf
            for aⱼ ∈ atomlist
                # 原子 aᵢ へ 結合する
                mol₂ = add_atom(mol, aⱼ)
                # 原子indexは最後
                j = natom(mol₂)
                bond = Bond(1)
                add_bond!(mol₂, i, j, bond)
                push!(mols, mol₂)
            end
        end
    end
    return mols
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
        for order₂ ∈ bonds
            bond₂ = Bond(order₂)
            push!(mols, add_bond(mol, i, j, bond₂))
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
        for orderᵢⱼ ∈ bonddict[order]
            bondᵢⱼ = Bond(orderᵢⱼ)
            if orderᵢⱼ == 0
                # 結合削除時 1  -> 0のときは
                # edge削除
                mol₂ = remove_bond(mol, i, j)
                push!(mols, mol₂)
            else
                push!(mols, add_bond(mol, i, j, bondᵢⱼ))
            end
        end
    end
    return mols
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
    addable_atom(mol)
end

function example3()
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
    addable_bond(mol)
    removable_bond(mol)
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

function example4()
    mol1, mol2 = examplemol1()
    @show tomatrix(mol1)
    @show tomatrix(mol2)
end

function example5()
    mol1, mol2 = examplemol1()
    return mol1 == mol2
end

function examplemol1()
    mol = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond1 = Bond(2)
    bond2 = Bond(1)
    add_bond!(mol, 1, 2, bond1)
    add_bond!(mol, 2, 3, bond2)

    mol2 = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol2, C)
    add_atom!(mol2, C)
    add_atom!(mol2, C)
    bond1 = Bond(1)
    bond2 = Bond(2)
    add_bond!(mol2, 1, 2, bond1)
    add_bond!(mol2, 2, 3, bond2)

    return mol, mol2
end

function examplemol2()
    mol = Molecule()
    C = Atom(:C, 6, 4)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    add_atom!(mol, C)
    bond1 = Bond(1)
    bond2 = Bond(1)
    bond3 = Bond(1)
    bond4 = Bond(1)
    bond5 = Bond(1)
    bond6 = Bond(1)
    bond7 = Bond(2)
    # 1 loop 1 分岐
    add_bond!(mol, 1, 2, bond1)
    add_bond!(mol, 2, 3, bond2)
    add_bond!(mol, 3, 4, bond3)
    add_bond!(mol, 4, 5, bond4)
    add_bond!(mol, 5, 6, bond5)
    add_bond!(mol, 6, 1, bond6)
    add_bond!(mol, 2, 7, bond7)
    return mol
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
    # 注目している原子
    current_atom_index = -1
    current_bond_order = 1
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
            # 芳香環は別関数で対処する
        end
        # 原子と一致
        if haskey(atomdict, t)
            println("原子")
            # 原子を追加したら注目している原子は追加した原子
            @show current_atom_index, current_bond_order
            current_atom_index = natom(mol)
        elseif haskey(bonddict, t)
            println("結合")
            # 次に追加する原子の結合字数
            current_bond_order = bonddict[t]
        end
        next = iterate(x, i)
    end
    return mol, nothing
end

"""
Molecule型へ変換する
"""
function Base.convert(::Type{Molecule}, m::GraphMol)
    nodes = m.nodeattrs
    edges = m.edges
    edgeattrs = m.edgeattrs

    mol = Molecule()
    # 原子を追加
    for node ∈ nodes
        index = findfirst(a->a.name==node.symbol, atomlist)
        atom = atomlist[index]
        add_atom!(mol, atom)
    end
    # 結合を追加
    for (edge, edgeattr) ∈ zip(edges, edgeattrs)
        i, j = edge
        bond = Bond(edgeattr.order)
        add_bond!(mol, i, j, bond)
    end
    return mol
end

"""
GraphMol型へ変換する
"""
function Base.convert(::Type{GraphMol}, m::Molecule)
    # graphmol(edges, nodes, edgeattrs)
    A = atoms(m)
    A₂ = map(A) do atom
        isaromatic = false
        SmilesAtom(atom.name, 0, 1, nothing, isaromatic, :unspecified)
    end
    # SmilesAtom(:C, 0, 1, nothing, false, :unspecified)
    B, bondmaps = bonds(m)
    bonds₂ = map(B) do bond
        SmilesBond(bond.order, false, :unspecified, :unspecified)
    end
    es₂ = map(bondmaps) do bondmap
        (bondmap.src, bondmap.dst)
    end
    graphmol(es₂, A₂, bonds₂)
end

"""
    Base.:(==)(mol₁::Molecule, mol₂::Molecule)

分子の等価
# 比較する
- bond order
"""
function Base.:(==)(mol₁::Molecule, mol₂::Molecule)
    setss = permuteindex(mol₁)
    sets₂ = tomatrix(mol₂)
    isequal = any(map(sets₁ -> sets₁ == sets₂, setss))
    return isequal
end

"""
    tomatrix

MolecularのSet表現
各原子について以下のリストを求め，Setにする
原子名, 原子index i, 隣接原子index j, 結合次数ij
"""
function tomatrix(m::Molecule)
    molset = Set{Set}()
    n = natom(m)
    for i ∈ 1:n
        aᵢ = get_atom(m, i)
        indices = neighbors(m.graph, i)
        subset = Set()
        for j ∈ indices
            oᵢⱼ = bondorder(m, i, j)
            data = [aᵢ.name, i, j, oᵢⱼ]
            push!(subset, data)
            # @show aᵢ.name, i, j, oᵢⱼ
        end
        push!(molset, subset)
    end
    return molset
end

"""
    permuteindex(m::Molecule)

グラフのvertex indexを入れ替えた全セットを返す
グラフの等価性を求めるために使用する
"""
function permuteindex(m::Molecule)
    # indexいれかえ
    n = natom(m)
    # [1, n]の順列
    # n!通り
    setss = []
    sets = tomatrix(m)
    for indices ∈ permutations(1:n)
        sets₂ = deepcopy(sets)
        # 単射
        for sᵢ ∈ sets
            for sᵢⱼ ∈ sᵢ
                sᵢⱼ[2] = indices[sᵢⱼ[2]]
                sᵢⱼ[3] = indices[sᵢⱼ[3]]
            end
        end
        push!(setss, sets₂)
    end
    return setss
end

"""
    moltosmiles(m::Molecule)

Molecule -> SMILES
"""
function moltosmiles(m::Molecule)
    # order = 2 -> "="
    # order = 3 -> "#"
    # 結合が複数ある場合は(で分岐
    m = deepcopy(m)
    # isparsed
    for i ∈ 1:natom(m)
        set_prop!(m.graph, i, :isparsed, false)
    end
    # loop indicesをセットする
    m = closeloop(m)
    moltosmiles(m, -1, 1)
end

hasnext(indices, state) = iterate(indices, state) !== nothing

"""
`m`は分子
`i₀`は一つ前の対象原子
`i`は対象の原子
"""
function moltosmiles(m::Molecule, i₀::Int, i::Int)
    set_prop!(m.graph, i, :isparsed, true)
    aᵢ = get_atom(m, i)
    print(aᵢ.name)
    indices = neighbors(m.graph, i)
    lidxsᵢ = has_prop(m.graph, i, :loopidxs) ? Set(get_prop(m.graph, i, :loopidxs)) : Set()
    print(join(string.(lidxsᵢ), ""))
    indices = filter(indices) do x
        if has_prop(m.graph, x, :loopidxs)
            lidxsⱼ = Set(get_prop(m.graph, x, :loopidxs))
            # 共通idxがある場合は省く
            if (0 < length(lidxsᵢ ∩ lidxsⱼ))
                return false
            end
        end
        # loopidxsが同じのが含まれる場合はスルー
        return x ≠ i₀
    end
    if indices === nothing || length(indices) == 0
        return nothing
    end
    # 隣接原子数
    next = iterate(indices)
    index, state = next
    while next !== nothing
        j, state = next
        # loopindexがあったらつける
        if has_prop(m.graph, j, :loopidxs)
            lidxs = get_prop(m.graph, j, :loopidxs)
        end
        # isparsedな原子の場合は終了
        isloop = get_prop(m.graph, j, :isparsed)
        # @show state, i, j, isloop
        if isloop
            next = iterate(indices, state)
            continue
        end
        # 最後の結合ならば，()なしで表示
        if !hasnext(indices, state)
            printsmiles(m, i, j)
        else
            print("(")
            printsmiles(m, i, j)
            print(")")
        end
        next = iterate(indices, state)
    end
end

function printsmiles(m::Molecule, i::Int, j::Int)
    order = bondorder(m, i, j)
    if order == 2
        print("=")
    elseif order == 3
        print("#")
    end
    aⱼ = get_atom(m, j)
    # 再帰呼び出し
    moltosmiles(m, i, j)
end

function example6()
    m1 = examplemol2()
    moltosmiles(m1)
end

function example7()
    m1 = examplemol2()
    # moltosmiles(m1)
    println()
    closeloop(m1)
end

"""
graph vertex propの
vectorへpushする
"""
function push_prop!(g, i::Int, key::Symbol, value::Any)
    v = has_prop(g, i, key) ? get_prop(g, i, key) : []
    push!(v, value)
    set_prop!(g, i, key, v)
end

"""
最小閉路探索
vertexに情報を埋め込む
"""
function closeloop(m::Molecule)
    closelist = []
    dfs(m, 1, [], [], closelist)
    closelist2 = []
    for c ∈ closelist
        # 同一閉路は省く
        if c ∉ closelist2 && reverse(c) ∉ closelist2
            push!(closelist2, c)
        end
    end
    sort!(closelist, by = c -> length(c))
    for cᵢ ∈ closelist
        sᵢ = Set(cᵢ)
        for cⱼ ∈ closelist2
            sⱼ = Set(cⱼ)
            if sⱼ ≠ sⱼ ∩ sᵢ
                push!(closelist2, cᵢ)
                break
            end
        end
    end
    # ループインデックスを埋め込む
    m₂ = deepcopy(m)
    for (i, c) ∈ enumerate(closelist2)
        push_prop!(m₂.graph, c[begin], :loopidxs, i)
        push_prop!(m₂.graph, c[end-1], :loopidxs, i)
    end
    return m₂
end

"""
Deep-first search
find loops
"""
function dfs(m::Molecule, v::Int, isvisited, target, closelist)
    push!(isvisited, v)
    target = copy(target)
    if v ∈ target
        push!(target, v)
        push!(closelist, target)
        return target
    end
    push!(target, v)
    indices = neighbors(m.graph, v)
    for i ∈ indices
        # 直前の頂点以外
        i₀ = 1 < length(target) ? target[end-1] : -1
        if i ≠ i₀
            dfs(m, i, isvisited, target, closelist)
        end
    end
    return target
end