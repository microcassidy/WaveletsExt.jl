function dyadic_length(i::Int64)
    J = floor(log2(i))
    if exp2(J) != i
        throw(ArgumentError("input needs to be a power of two."))
    else
        J
    end
end
