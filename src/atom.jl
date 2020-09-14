export Atom

"""
`number`は原子番号で`valence`は原子価数
"""
struct Atom
    name::Symbol
    number::Int
    valence::Int
end

"""
原子価数
"""
function nvalence(a::Atom)::Int
    return a.valence
end

atomdata = [
    [:H, 1, 1],
    [:B, 5, 3],
    [:C, 4, 4],
    [:N, 7, 3],
    [:O, 8, 2],
    [:F, 9, 1],
]
atomlist = map(l->Atom(l...), TransMol.atomdata)