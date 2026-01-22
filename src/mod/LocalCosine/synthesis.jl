function synthesis_step!(out::AbstractVector{Float64},
                        v::AbstractVector{Float64},
                        lcb::LocalCosineBasis,
                        block_size::Int64,j::Int64)
end

function synthesis_operator(lcb::LocalCosineBasis,coef_tree::AbstractVector{Float64})
    N,_,nblocks = size(coef_tree)
    out = lc{Float64}(undef,N*nblocks)
    synthesis_operator!(out,lcb::LocalCosineBasis,coef_tree::Array{Float64,3})
end

function synthesis_operator!(out::Vector{Float64}, coef_tree::Array{Float64,3},
                            lcb::LocalCosineBasis,
                            v::AbstractVector{Float64})
    #TODO: add logic for best basis
    N,J,nblocks = size(coef_tree)
    for blockidx in 1:nblocks
        for j in 0:J-1
            vhat = view(coef_tree,:,j+1,block_idx)
            error("unimplemented")
            synthesis_step!(out)

        end
    end
end
