#bp is interior of window
#bm is the [m]irror
function unfoldedge!(::Type{LeftPacket}, lcb::LocalCosineBasis)
	m = lcb.m
	front = 1:m
	# for idx in 0:lcb.m-1
		# lcb.packet.left[begin+idx] = lcb.packet.centre[begin+idx] * (1 - lcb.bell.interior[begin+idx]) / lcb.bell.interior[begin+idx]
	lcb.packet.left[front] = lcb.packet.centre[front] .* (1 .- lcb.bell.interior) ./ lcb.bell.interior
	# end
end
function unfoldedge!(::Type{RightPacket}, lcb::LocalCosineBasis)
	# set!(LeftPacket, lcb, lcb.centre)
	m = lcb.m
    @info "unfold-lasindex centre : $(lastindex(lcb.packet.centre))"
	n = length(lcb.packet.centre)
	@info "n: $n"
	# back  = n:-1:(n-m+1)
	back  = lastindex(lcb.packet.centre):-1:lastindex(lcb.packet.centre)-m+1
	front = 1:m
	# lcb.packet.centre[end-m+1:end] += lcb.packet.centre[end-idx]

	lcb.packet.right[front] = lcb.packet.centre[back] .* (1 .- lcb.bell.interior) ./ lcb.bell.interior

	# lcb.packet.centre[]+= lcb.packet.centre[back] .* (1 .- lcb.bell.interior) ./ lcb.bell.interior
	# nothing
end
	# for idx in 0:m-1
	# 	# lcb.packet.right[end-idx] *= lcb.packet.centre[end-idx]*(1 - lcb.bell.interior[idx+begin]) / lcb.bell.interior[idx+begin]
	# 	lcb.packet.centre[end-m+1+idx] +=
	# end
# end
	# if strcmp(which,'left'),
	# 	front = 1:m;
	# 	extra(front) = xc(front) .* (1-bp)./bp;
	# else
	# 	back  = n:-1:(n+1-m);
	# 	extra(back) = xc(back) .* (1-bp)./bp;
	# end	

function unfold!(lcb::LocalCosineBasis)
	m = lcb.m
	n = length(lcb.packet.centre)
	lcb.packet.right .= 0
	lcb.packet.left .= 0
	back  = n:-1:(n-m+1)
	front = 1:m

	lcb.packet.left = lcb.bell.exterior[reverse(front)] .* lcb.packet.centre[front]
	# lcb.packet.right= (-1) .* reverse(lcb.bell.exterior) .* lcb.packet.centre[back]
	lcb.packet.right= (-1) .* lcb.bell.exterior[reverse(front)] .* lcb.packet.centre[back]

	lcb.packet.centre[front] .*=  lcb.bell.interior
	lcb.packet.centre[back]  .*=  lcb.bell.interior

#  for i in 0:lcb.m-1
# 		lcb.packet.left[end-i] = lcb.packet.centre[i+begin] * lcb.bell.exterior[end-i]
# 		lcb.packet.right[i+begin] = -lcb.packet.centre[end-i] * lcb.bell.exterior[end-i]
# #UNCSURE OF ABOVE
# 		lcb.packet.centre[i+begin] *= lcb.bell.interior[begin+i]
# 		lcb.packet.centre[end-lcb.m+1+i] *= lcb.bell.interior[end-i]
# 	end
end
	# n = length(y);
	# m = length(bp);
	# xc = y;
	# xl = 0 .*y;
	# xr = 0 .*y;
	# front = 1:m;
	# back  = n:-1:(n+1-m);
	# xc(front) =       bp .* y(front);
	# xc(back)  =       bp .* y(back );
	# xl(back)  =       bm .* y(front); <<<<<<<<<<<<

	# xr(front) =     (-1) .* bm .* y(back);


# xl = bm .* reverse(y)

# reverse(left) = reverse(bm) .* yf
# left = bm * reverse(y)
