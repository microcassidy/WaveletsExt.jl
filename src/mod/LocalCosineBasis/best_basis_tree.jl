
function basis_cost(tr,ma)
    costs = Vector{Float64}(undef, sum(2^j for j in 0:max_depth))
    i = 1
    for j in 0:max_depth
        for xs in blocks
            for x in xs
                cost[i] += x^2
            end
            i+=1
        end
    end
    @assert i == length(costs) + 1
end
