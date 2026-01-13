

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
function _fold!(h::VariableVector,centre::AbstractVector{Float64},
                  left::AbstractVector{Float64},
                  right::AbstractVector{Float64},
                  bell::OrthonormalBell)
    m = length(bell)
    n = length(centre)

    @info "lengths: bell $m,centre $n,left: $(length(left)), right $(length(right))"
    @assert length(h) == length(centre)

    # out = similar(centre)
    # fill!(out,0)
    @info size(centre)
    for idx in 0:m-1
        #front
        h[idx+begin]= interior(bell)[idx+begin] * centre[idx+begin] + exterior(bell)[idx+begin] * left[end-m+1+idx]
        #back
        @debug "h idx $(length(h) - m + 1 + idx)"
        h[end-m+1+idx]  = interior(bell)[begin+idx] * centre[end-m+1+idx] - exterior(bell)[begin+idx] * right[begin+idx]
    end
    nothing
end
