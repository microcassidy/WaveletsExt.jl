
# function synthesis_operator(lcb::LocalCosineBasis,v::Vector{Float64},basis_tree::BitVector)
#     N,_,nblocks = size(coef_tree)
#     out = lc{Float64}(undef,N*nblocks)
#     synthesis_operator!(out,lcb::LocalCosineBasis,coef_tree::Array{Float64,3})
# end

export SynthesisPlan
import .Utils:BinaryTree
struct SynthesisPlan <: AbstractVector{Tuple{UnitRange,Int64}}

    subsets::Vector{UnitRange}
    dctidx::Vector{Int64}
    function SynthesisPlan(tr::BitVector,N::Int64)
        idxs = findall(!iszero,tr)
        subsets = similar(idxs,UnitRange)
        dctidx = map!(x->log2(x)|>floor|>Int64,idxs) #also are the js
        from = 1
        for idx in eachindex(dctidx)
            sz = div(N,2^dctidx[idx])
            subsets[idx] = from:sz-1
            from += sz
        end
        new(subsets,dctidx)
    end
end
Base.getindex(s::SynthesisPlan,i::Int)::Tuple{UnitRange,dctidx} = (s.subsets[i],s.dctidx[i])
Base.length(s::SynthesisPlan,i::Int) = length(s.subsets)
Base.size(s::SynthesisPlan) = length(s.subsets)
Base.firstindex(s::SynthesisPlan) = 1
Base.lastindex(s::SynthesisPlan) = length(s.subsets)


function synthesis_step!(out::AbstractVector{Float64}, v::AbstractVector{Float64}, lcb::LocalCosineBasis, synthesis_plan::SynthesisPlan)
    for idx in eachindex(synthesis_plan)
        unitrng,dct_idx = synthesis_plan[idx]
        F = @view lcb.dct_plans[dct_idx]
        lcb.packet.centre = F \ v[unitrng]
        # idx == firstindex(idx) ? unfoldedge!(LeftPacket,lcb,v[synthesis_plan[idx][1]]) ? set!(LeftPacket,v[synthesis_plan[idx-1][1]])
        # idx == firstindex(idx) ? unfoldedge!(LeftPacket,lcb,v[synthesis_plan[idx][1]]) ? set!(LeftPacket,v[synthesis_plan[idx-1][1]])

    end
end


pairs(s::SynthesisPlan) = zip(s.subsets, s.dctidx)


lsum!(p::Packet) = p.centre[1:length(p.left)] += p.left
rsum!(p::Packet) = p.centre[1:length(p.right)] += p.right

lsum!(b::LocalCosineBasis) = b.packet.centre[1:length(b.packet.left)] += b.packet.left
rsum!(b::LocalCosineBasis) = b.packet.centre[end-b.m+1:end] += b.packet.right
function lrsum(p::Packet)
    lsum(p)
    rsum(p)
end

function synthesis_operator!(out::Vector{Float64}, v::Vector{Float64},
                            lcb::LocalCosineBasis,plan::SynthesisPlan,bb::BitVector)
    # N,J,nblocks = size(c)
    # N =
    # isleft = true
    # offset=0
    # subsets =
    basis_idxs = findall(bb)
    nblocks_bb = length(bb) #TODO:unsure

    nblocks_sig = Int64(length(v) / lcb.N)


    col_idxs = map(idx->Int64(floor(log2(idx)))+1,basis_idxs)
    @info col_idxs
    row_ranges = map(idx->getrowrange(BinaryTree,lcb.N,idx),basis_idxs)
    # col_ranges = map(idx->getcolrange(lcb.N,idx,:binary),basis_idxs)
    for block_offset in 0:nblocks_sig-1
            for (count,(ri,ci)) in enumerate(zip(row_ranges,col_idxs))
                # @assert ri
                vw = @view v[block_offset*lcb.N+begin:(block_offset+1)*lcb.N]
                ov = @view out[block_offset*lcb.N+begin:(block_offset+1)*lcb.N]
                fill!(lcb.packet.left,0)
                fill!(lcb.packet.right,0)
                _update_size!(lcb.packet,length(ri))

                F = lcb.dct_plans[ci]
                @info F
                @info length(ri)
                set!(CentrePacket,lcb,F \ vw[ri])
                unfold!(lcb)
                count == 1 ? unfoldedge!(LeftPacket,lcb) : lsum!(lcb)
                count == length(basis_idxs) ? unfoldedge!(RightPacket,lcb) : rsum!(lcb)
                ov[ri] = lcb.packet.centre
            end
    end


        # fi,li = first(vrange), last(vrange)
        # out[fi-m+1:fi-1] += idx==firstindex(plan) ? unfoldedge!(LeftPacket,lcb) : lcb.packet.left
        # out[li+1:li+m] += idx==firstindex(plan) ? unfoldedge!(RightPacket,lcb) : lcb.packet.right
        # idx==lastindex(plan) && unfoldedge!(RightPacket,v[vrange])
        # copyto!(out[rng],packet.centre[begin:end])
end
