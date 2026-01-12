
mutable struct Folded <: AbstractVector{Float64}
    data::Vector{Float64}
    length::Int
    max_length::Int
    Folded(m::Int) = new(zeros(m),m,m)
end
_update_size!(v::Folded,m::Int) = 1<=m<=v.max_length ? setfield(v,:length,m) : error("outside bounds of 1 and $(v.max_length)")
Base.length(v::Folded) = v.length
Base.size(v::Folded) = v.length
Base.getindex(v::Folded,i::Int) = v.data[i]

function pseudopacket(centre_packet::AbstractVector{Float64},bell::OrthonormalBell, side::Symbol)
    """
    produces a pseudopacket so that the folding function doesn't do anything
    TODO: check whether this function could just be a noop
    """
    packet = zeros(length(centre_packet))
    N = length(packet)
    m = length(bell)
    if side == :left
        for idx in 0:m-1
            packet[end - N + 1 + idx] = centre_packet[idx+begin] * (1 - interior(bell)[idx+begin]) / exterior(bell)[idx+begin]
        end
    elseif side == :right
        for idx in 0:m-1
            packet[idx+begin] = -centre_packet[end-N+1+idx] * ( 1 - interior(bell)[idx+begin]) / exterior(bell)[idx+begin]
        end
    else
        error("side should be $(:left) or $(:right)")
    end
    packet
end
function _fold!(h::Folded,centre::AbstractVector{Float64},
                  left::AbstractVector{Float64},
                  right::AbstractVector{Float64},
                  bell::OrthonormalBell)
    m = length(bell)
    n = length(centre)

    # out = similar(centre)
    # fill!(out,0)
    for idx in 0:m-1
        #front
        h[idx+begin]= interior(bell)[idx+begin] * centre[idx+begin] + exterior(bell)[idx+begin] * left[end-m+1+idx]
        #back
        h[end-m+1+idx]  = interior(bell)[idx+begin] * centre[end-m+1+idx] - exterior(bell)[idx+begin] * right[idx+begin]
    end
    out
end
