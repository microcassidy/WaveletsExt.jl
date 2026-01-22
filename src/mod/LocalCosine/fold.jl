function foldedge!(::Type{LeftPacket}, lcb::LocalCosineBasis)
    """
    Notes:
    left: f(-t) = f(t)(1 - g(t)) / g(-t) on [0,a]
        - The outer reversed to account for the filter direction evaulation (during summation)
        - the inner filter needs reversing for the same reason
    """
    m = lcb.m

    front = view(lcb.packet.centre,1:m)

    # the outer reverse is for correct assignment, the innner is for correct filter direction
    # The filter is specified on the left only
    lcb.packet.left .*= lcb.packet.centre[1:m]
    t = reverse((front .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior))) #unsimplified
    @assert length(t) == lcb.m
    set!(LeftPacket,lcb, reverse(front) .* (1 .- reverse(lcb.bell.interior)) ./ (lcb.bell.exterior))
    @assert t == lcb.packet.left "something went wrong with left packet assumption"
    nothing
end
function foldedge!(::Type{RightPacket}, lcb::LocalCosineBasis)
    """
    right: f(2-t) = -f(t)/g(2-t) * (1 - g(t)) on [1-a, 1]
         - the interior filter needs reversing as it is for the LHS
         - The time reversed exterior filter is the forward filter rev(rev(x)) = x and therefore does not need to be flipped
         - The time reversed exterior filter is the forward filter rev(rev(x)) = x and therefore does not need to be flipped
         - the above also makes sense because the sum is contracting  >---1---< as we are moving from left to right
            - it is accounting for the imaginary rising edge of it's neighbour

         - The entire statement needs reversing as it is for the values of the time reversed component on the exterior
         - If the outer reverse is 'evaluated' it will simplifify to reverse(back) .* (1 .- b.interior) ./ reverse(b.exterior)??
            -an imaginary rising edge of the exterior
    """
    m = lcb.m
    back = view(lcb.packet.centre,lastindex(lcb.packet.centre)-m+1:lastindex(lcb.packet.centre))
    t = - reverse(back .* (1 .- reverse(lcb.bell.interior)) ./ lcb.bell.exterior) #unsimplified statement
    set!(RightPacket,lcb, -reverse(back) .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior)) #simplified)
    @assert t == lcb.packet.right "something went wrong with right packet assumption"
    # copyto!(p.right, -reverse(back .* (1 .- U.interior) ./ reverse(U.exterior))) #TODO:check
end


function fold!(lcb::LocalCosineBasis)
    m = lcb.m
    # @inline front_m(v::AbstractVector{float}) =
    # @inline back_m(v::AbstractVector{float}) =
    # front = view(lcb.packet.centre,firstindex(lcb.packet.centre):m)
    # back = view(lcb.packet.centre,lastindex(lcb.packet.centre)-m+1:lastindex(lcb.packet.centre))

    #left edge
    # h(t) = f(t)g(t) + f(-t)g(-t) on [0,a]
    for i in 1:m
        lcb.packet.centre[i] *= lcb.bell.interior[i]
        lcb.packet.centre[i] += lcb.bell.exterior[end-lcb.m+i] * lcb.packet.right[end-lcb.m+i]
        lcb.packet.centre[end-lcb.m+i] *= lcb.bell.interior[end-lcb.m+i]
        lcb.packet.centre[end-lcb.m+i] -= lcb.bell.exterior[i] * lcb.packet.right[end-lcb.m+i] #reversed as the filter is on the right, the filter doesn't need reversing twice
    end
    # front .*= lcb.bell.interior
    # front += reverse(lcb.bell.exterior .* lcb.packet.right) #needs to be reversed so the terms are summed <---outwards---->
    #right edge
    # #h(t) = f(t)g(t) - f(2-t)g(2-t)
    # back *= reverse(lcb.bell.interior) #reversed as the filter mirrored on the right
    # #The inner reverse is for the mirror filter and the outer is for time reversal
    # back -= reverse(reverse(U.exterior) .* packet.right) #reversed as the filter mirrored on the right
    # back -= lcb.bell.exterior .* reverse(lcb.packet.right) #reversed as the filter is on the right, the filter doesn't need reversing twice

    nothing
end
