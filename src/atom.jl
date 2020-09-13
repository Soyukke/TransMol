export Atom

"""
`number`は原子番号で`valence`は原子価数
"""
struct Atom
    number::Int
    valence::Int
end

"""
原子価数
"""
function nvalence(a::Atom)::Int
    return a.valence
end