function cost(coef_tree::Array{Float64,3})
    n,m,_ = size(coef_tree)
    out = zeros(n,m)
    cost!(out,coef_tree)
end

function grouped_ssq!(o::AbstractVector{Float64}, v::Vector{Float64},j::Int64)
    n = length(v)
    block_size = div(n,2^j)
    nblocks = 2^j
    bs = reshape(v,block_size,nblocks)
    for idx in 1:nblocks
        o[idx] += norm(bs[:,idx],2)^2
    end
end

export cost!
function cost!(out::Vector{Float64},m::AbstractMatrix{Float64})
    N,J = size(m)
    bottom = 1
    for j in 0:J-1
        block_size = div(N,2^j)
        nblocks = 2^j

        ov = @view out[bottom:bottom+nblocks-1]
        grouped_ssq!(ov,m[:,j+begin],j)
        bottom += nblocks
    end
end

"""
expectation of the sum of squares
recursive calculcation ``\\mu_{n+1} = \\frac{n}{n+1}\\mu_{n} + \\frac{x_n}{n+1}``
"""
function cost!(out::Vector{Float64}, coef_tree::Array{Float64,3})
    _,_,N = size(coef_tree)
    T = eltype(out)
    @views for n in 0:N-1
        out .*= T(n)
        cost!(out, coef_tree[:,:,n+begin])
        out ./= T(n+1)
    end
end

export cost_test
function cost_test(coef_tree)
    nlevels = size(coef_tree,2)-1
    nblocks = size(coef_tree,3)
    out = zeros(2^(nlevels+1) - 1)
    for x in eachslice(coef_tree,dims=3)
        ss = norm(x[:,1],2) #total energy
        s
        bottom = 1
        for (idx,c) in enumerate(eachcol(x))
            m = reshape(c,:,2^(idx-1))
            t  = abs.(x) / ss #normalize
            v = vcat(map(x->sum(t.^2),eachcol(m)))
            @info "vlength :$(length(v))"
            out[bottom:bottom+length(v)-1] += v
            bottom+=length(v)
        end
    end
    return out ./ size(coef_tree,3)
end
