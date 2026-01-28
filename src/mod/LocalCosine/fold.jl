"""
Notes:
left: f(-t) = f(t)(1 - g(t)) / g(-t) on [0,a]
    - The outer reversed to account for the filter direction evaulation (during summation)
    - the inner filter needs reversing for the same reason
"""
function foldedge!(::Type{LeftPacket}, lcb::LocalCosineBasis)
    m = lcb.m
	n = length(lcb.packet.centre);
	back  = n:-1:(n-m+1)
	front = 1:m

    lcb.packet.left[reverse(front)] = lcb.packet.centre[front] .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior)
    nothing
end

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
function foldedge!(::Type{RightPacket}, lcb::LocalCosineBasis)
    m = lcb.m
	n = length(lcb.packet.centre);
	back  = n:-1:(n-m+1)
	front = 1:m
    #
    lcb.packet.right = (-1) .* lcb.packet.centre[back] .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior)
    nothing
    # # t = - reverse(back .* (1 .- reverse(lcb.bell.interior)) ./ lcb.bell.exterior) #unsimplified statement
    # # # set!(LeftPacket,lcb, -reverse(back) .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior)) #simplified)
    # # set!(RightPacket,lcb, -reverse(back) .* (1 .- lcb.bell.interior) ./ reverse(lcb.bell.exterior)) #simplified)
    # # lcb.packet.right = -lcb.packet.centre[back] .* (1 .- lcb.bell.interior)./reverse(lcb.bell.exterior)

    # for idx in 0:m-1
    #     # lcb.packet.right[end-idx] = -lcb.packet.centre[end-m+1+idx]
    #     # lcb.packet.right[end-idx] *= (1 - lcb.bell.interior[end-idx])
    #     # lcb.packet.right[end-idx] /= lcb.bell.exterior[begin+idx]

    #     # lcb.packet.right[end-idx] = -lcb.packet.centre[end-m+1+idx]*(1 - lcb.bell.interior[end-idx])/lcb.bell.exterior[begin+idx]
    #     lcb.packet.right[end-idx] = 0
    # end

    # @assert t == lcb.packet.right "something went wrong with right packet assumption"
    # copyto!(p.right, -reverse(back .* (1 .- U.interior) ./ reverse(U.exterior))) #TODO:check
end


function fold!(lcb::LocalCosineBasis)
    m = lcb.m
    for i in 0:m-1
        lcb.packet.centre[begin+i] *= lcb.bell.interior[begin+i]
        lcb.packet.centre[begin+i] += lcb.bell.exterior[end-i] * lcb.packet.left[end-i]

        lcb.packet.centre[end-m+1+i] *= lcb.bell.interior[end-i]
        lcb.packet.centre[end-m+1+i] -= lcb.bell.exterior[begin+i] * lcb.packet.right[end-i] #reversed as the filter is on the right, the filter doesn't need reversing twice
    end

    nothing
end
