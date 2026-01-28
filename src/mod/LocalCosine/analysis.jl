const DEFAULT_BASIS_TOPLEVEL_LENGTH = 1024
fftwsimilar(sz) = Array{Float64}(undef,sz...)


function plan_dct_iv(n,max_depth::Integer)::Vector{r2rFFTWPlan}
    plans::Vector{DCT_T} = [plan_r2r(fftwsimilar(n ÷ 2^j),REDFT11) for j in 0:max_depth]
end

function coefficient_tree(n::Int64,max_depth::Int64)
    Matrix{Float64}(undef, max_depth + 1)
end
# reset!(tr::Matrix{Float64}) = fill!(tr,0)

function dyadic_length(i::Int64)
    J = floor(log2(i))
    exp2(J) == i ? J : error("vector length needs to be a power of two.")
end



# import Base: *
# *(F::FFTWPlan,lcb::LocalCosineBasis) = F * lcb.packet.centre

fix_n(f,args...) = foldl(Base.Fix1,args;init=f)
"""
out and packet are modified, the first 3 in place for
partial application convenience
"""
# dctiv!(out::AbstractVector{})

function analysis_step!(out::AbstractVector{Float64},
                        v::AbstractVector{Float64},
                        lcb::LocalCosineBasis,
                        block_size::Int64,j::Int64)
    endi = 2^j-1
    f = lcb.dct_plans[j+1]
    # @boundscheck checkbounds(v,1:block_size*endi)
    # @info "packet length $(length(lcb.packet.centre))"
    m = lcb.m
    ranges = [block_size*i+1:block_size*(i+1) for i in 0:endi]
    for idx in 0:endi
        crange = ranges[idx+begin]
        rrange=nothing
        lrange=nothing
        copyto!(lcb.packet.centre,v[crange]);
        if idx == 0
            foldedge!(LeftPacket, lcb)
        else
            lrange = ranges[idx-1+begin]
            lrange = lrange.stop-m+1:lrange.stop
            lcb.packet.left = v[lrange]
        end

            # set!(LeftPacket,lcb,v[(idx-1)*block_size+begin:(idx)*block_size])
        if idx == endi
            foldedge!(RightPacket,lcb)
        else
            # trange = idx*block_size+firstindex(v):(idx+1)*block_size
            rrange = ranges[idx+1+begin]
            rrange = rrange.start:rrange.start+m-1
            # #@info "$lrange $crange $rrange"
            lcb.packet.right =v[rrange]
            # set!(RightPacket,lcb)
        end
        fold!(lcb)
        out[(idx*block_size)+1:(idx+1)*block_size] = f * lcb.packet.centre
    end
end

function analysis_operator(lcb::LocalCosineBasis,v::AbstractVector{Float64})
    nblocks = Int64(div(length(v),lcb.N))
    # nblocks * lcb.N != lcb.N && @warn "v is not a multiple of the block size nblocks=$nblocks $(nblocks * lcb.N)"
    # nblocks * lcb.N != lcb.N
    coef_tree = Array{Float64}(undef, (lcb.N,lcb.J + 1, nblocks)) #(n, n_decompositions, n_blocks)
    analysis_operator!(coef_tree,lcb,v)
end
function analysis_operator!(coef_tree::Array{Float64,3}, lcb::LocalCosineBasis, v::AbstractVector{Float64})
    n,m,o = size(coef_tree)
    #@info size(coef_tree)
    # @boundscheck checkbounds(v, 1:n*o,v)
    @assert m - 1 == lcb.max_depth "$m:$(lcb.max_depth)"
    for blockindex in 0:o-1
        for j in 0:lcb.max_depth
            reset!(lcb.packet) #probably unnecessary, just as a precaution
            block_size = lcb.N ÷ 2^j
            _update_size!(lcb.packet,block_size)
            analysis_step!(view(coef_tree,:,j+1, blockindex+1),
                           v[(blockindex)*lcb.N+begin:lcb.N*(blockindex + 1)],
                           lcb,
                           block_size, j)
            # #@info "---"
            # # for i in 0:j
            # #@info coef_tree[1:4,j+begin]
            # # end
            # #@info "---"
            # j == 2 && error("foo")
        end
    end
end
