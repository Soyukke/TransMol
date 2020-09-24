"""
SMILES Strings

https://www.daylight.com/cheminformatics/index.html
"""
smiles_tokens = [
    "B", "C", "N", "O", "S", "P", "F", "I", "Cl", "Br",
    "c", "n", "o", "s",
    "[", "]", "(", ")",
    "1", "2", "3", "4", "5", "6", "7", "8",
    "H", "-", "+",
    "=", "#",
    "@", "@@",
    "/", "\\", "."
]

atomtokens = ["B", "C", "N", "O", "S", "P", "F", "I", "Cl", "Br"]
aromatictokens = ["c", "n", "o", "s"]

"""
    aromatic2atom(token::String)

# Examples
```julia
julia> aromatic2atom("br")
"Br"
```
"""
function aromatic2atom(token::String)
    return uppercasefirst(token)
end

struct Smiles <: AbstractString
    smiles::String
    tokens::Vector{String}
end

function Smiles(str::String)
    strs = split(str, "")
    next = iterate(strs)
    tokens = String[]
    token = ""
    state = -1
    while next !== nothing
        # 2文字tokenの探索
        token₁, state₁ = next
        next₂ = iterate(strs, state₁)
        if next₂ !== nothing
            token₂, state₂ = next₂
            token₃ = token₁ * token₂
            if token₃ ∈ smiles_tokens
                token = token₃
                state = state₂
            else
                token = token₁
                state = state₁
            end
        else
            token = token₁
            state = state₁
        end
        if token ∈ smiles_tokens
            push!(tokens, token)
        else
            throw(UndefVarError(Symbol(token)))
        end
        # 1文字tokenの探索
        next = iterate(strs, state)
    end
    return Smiles(str, tokens)
end

"""

SMILESをトークン単位で取得する
"""
function Base.getindex(s::Smiles, i::Int)
    return s.tokens[i]
end

function Base.length(s::Smiles)
    return length(s.tokens)
end

function Base.show(io::IO, s::Smiles)
    println("Number of Tokens: ", length(s))
    print(io, "SMILES: ")
    Base.show(io, s.smiles)
end

Base.iterate(smiles::Smiles) = iterate(smiles, 1)
function Base.iterate(smiles::Smiles, i::Integer)
    if length(smiles) < i
        return nothing
    end
    return smiles[i], i + 1
end

"""
for debug
"""
function main()
    s = "C=CCBrO@@CCl"
    smi = Smiles(s)
    @show smi.tokens
    ""
    return smi
end