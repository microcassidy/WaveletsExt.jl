function cost(tr::Matrix{Float64})::Vector{Float64}
    nr,nc = size(tr)
    J = nc - 1
    ncost = (1 << nc) - 1
    cost_tr = zeros(ncost)
    i = 0
    for (j,c) in enumerate(eachcol(tr))
        nsplit = 2^(j-1)
        bsize = Int(nr / nsplit)
        @views for si in 0:nsplit-1
            v = c[si*bsize+1:(si+1)*bsize]
            cost_tr[i+begin]= v'v
            i+=1
        end
    end
    @assert i == ncost "ncost:$ncost i:$(i)"
    return cost_tr
end
