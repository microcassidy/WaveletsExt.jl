struct FoldingOperator
    length::Integer
    inner_support::Integer
    lskirt::Vector{Float64}
end

foldsimilar()= Array{T}
function *(F::FoldingOperator,v::AbstractVector{Float64})
    length(v) - length(lskirt) = inner_support
    mul!(out,F,v)
end
mirror(F::FoldingOperator) = view(F,length(F):1)
function mul!(dest,F,v)
    m = mirror(F)
    rskirt =
end
