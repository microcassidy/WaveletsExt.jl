
# function synthesis_operator(lcb::LocalCosineBasis,v::Vector{Float64},basis_tree::BitVector)
#     N,_,nblocks = size(coef_tree)
#     out = lc{Float64}(undef,N*nblocks)
#     synthesis_operator!(out,lcb::LocalCosineBasis,coef_tree::Array{Float64,3})
# end

export SynthesisPlan
import .Utils:BinaryTree
struct SynthesisPlan <: AbstractVector{Tuple{UnitRange,Int64}}
# import Datastructures:Stack

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
function rsum!(b::LocalCosineBasis)
    @info "lasindex centre : $(lastindex(b.packet.centre))"
    b.packet.centre[lastindex(b.packet.centre)-b.m+1:lastindex(b.packet.centre)] += b.packet.right
end
function lrsum!(b::LocalCosineBasis)
    lsum!(b)
    rsum!(b)
end


@inline _nodeidx(depth,block)::Int64 = (1 << depth) + block
function synthesis_operator!(out::Vector{Float64}, v::Vector{Float64},
                            lcb::LocalCosineBasis,bb::BitVector)
    # N,J,nblocks = size(c)
    # N =
    # isleft = true
    # offset=0
    # subsets =
    basis_idxs = findall(bb)
    nblocks_bb = length(bb) #TODO:unsure

    nblocks_sig = Int64(length(v) / lcb.N)


    col_idxs = map(idx->Int64(floor(log2(idx)))+1,basis_idxs)
    row_ranges = map(idx->getrowrange(BinaryTree,lcb.N,idx),basis_idxs)
    @assert length(col_idxs) == length(row_ranges)
    # col_ranges = map(idx->getcolrange(lcb.N,idx,:binary),basis_idxs)
    @debug "LOOP ENTRY-----"
    #TODO: the ordering is wrong and will break for more complex basis trees I think
    #need inorder traversal
    #
    # stack = Stack{Tuple{Int64,Int64}}([(0,0)])
    out .= 0


    @inline m = lcb.m
    for block_offset in 0:nblocks_sig-1
            for (count,(ri,ci)) in enumerate(zip(row_ranges,col_idxs))
                @debug "count: $count"
                @debug "\t--LOOP BEGIN---"
                @debug "ri:$ri ci:$ci"
                _update_size!(lcb.packet,length(ri))
                lcb.packet.left .=0
                lcb.packet.right .= 0
                lcb.packet.centre .= 0
                # @assert ri
                rng = block_offset*lcb.N+1:(block_offset+1)*lcb.N
                vw = @view v[rng]
                ov = @view out[rng]

                @info "ri length :$(length(ri))"
                @info (ri.start,ri.stop)

                @info (firstindex(vw),lastindex(vw))
                @info (firstindex(ov),lastindex(ov))

                fill!(lcb.packet.left,0)
                fill!(lcb.packet.right,0)
                @warn "in testing, do not use"
                F = lcb.dct_plans[ci]
                # set!(CentrePacket,lcb,F \ vw[ri])
                copyto!(lcb.packet.centre,F \ vw[ri])

                unfold!(lcb)

                ov[ri] += lcb.packet.centre

                if ri.start == 1
                    @debug "unfolding left packet"
                    unfoldedge!(LeftPacket,lcb)
                    # ov[ri.start:ri.start+m-1] += lcb.packet.left
                    ov[begin:m] += lcb.packet.left
                else
                    ov[reverse(ri.start-1-m+1:ri.start-1)] += lcb.packet.left
                    # lsum!(lcb)
                end

                # if count == lastindex(col_idxs)
                if ri.stop == lastindex(v)
                    @debug "unfolding right packet"
                    unfoldedge!(RightPacket,lcb)
                    ov[reverse(ri.stop-m+1:end)] += lcb.packet.right
                else
                ov[ri.stop+1:ri.stop+m] += lcb.packet.right
                # rsum!(lcb)
                end
                # lrsum!(lcb.packet)
                # ov[ri] = lcb.packet.centre
                @debug "\t--LOOP END---"
            end
    end


        # fi,li = first(vrange), last(vrange)
        # out[fi-m+1:fi-1] += idx==firstindex(plan) ? unfoldedge!(LeftPacket,lcb) : lcb.packet.left
        # out[li+1:li+m] += idx==firstindex(plan) ? unfoldedge!(RightPacket,lcb) : lcb.packet.right
        # idx==lastindex(plan) && unfoldedge!(RightPacket,v[vrange])
        # copyto!(out[rng],packet.centre[begin:end])
end
