export Atom

atomdata = [
    [:H, 1, 1],
    [:B, 5, 3],
    [:C, 4, 4],
    [:N, 7, 3],
    [:O, 8, 2],
    [:F, 9, 1],
    [:S, 16, 6],
]

"""
`number`は原子番号で`valence`は原子価数
"""
struct Atom
    name::Symbol
    number::Int
    valence::Int
end
atomlist = map(l->Atom(l...), TransMol.atomdata)

function Atom(name::Symbol)
    idx = findfirst(x -> x.name == name, atomlist)
    return atomlist[idx]
end

"""
原子価数
"""
function nvalence(a::Atom)::Int
    return a.valence
end

