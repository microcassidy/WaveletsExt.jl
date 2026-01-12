function dyad_length(v::Vector{<:Number})
    n = length(v)
    J = log2(n) |> Int
    return J,n
end
