function unfoldedge(::Type{LeftPacket}, lcb::LocalCosineBasis)
	# n = length(xc);
	# m = length(bp);
	# extra = xc.*0;

	for idx in 1:lcb.m
		lcb.packet.left[idx] = lcb.packet.centre[idx] * (1 - lcb.bell.interior[idx])/lcb.bell.interior[idx]
	end
end
function unfoldedge(::Type{RightPacket}, lcb::LocalCosineBasis)
	# n = length(xc);
	# m = length(bp);
	# extra = xc.*0;
	for idx in 0:lcb.m-1
		#TODO
		# lcb.packet.right[end - idx] = (1 - lcb.bell.exterior[e])lcb.bell.exterior
	end
	# if strcmp(which,'left'),
	# 	front = 1:m;
	# 	extra(front) = xc(front) .* (1-bp)./bp;
	# else
	# 	back  = n:-1:(n+1-m);
	# 	extra(back) = xc(back) .* (1-bp)./bp;
	# end
end
