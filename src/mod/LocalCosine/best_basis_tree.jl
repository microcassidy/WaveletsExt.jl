function cost(coef_tree::Array{Float64,3})
    n,m,_ = size(coef_tree)
    out = zeros(n,m)
    cost!(out,coef_tree)
end

function cost!(out::Matrix{Float64}, coef_tree::Array{Float64,3})
    """
    recursive calculcation $\mu_{n+1} = \frac{n}{n+1}\mu_{n} + \frac{x_n}{n+1}$
    """
    fill!(out,0)
    L,J,N = size(coef_tree)
    for n in 0:N-1
        for j in 1:J
            for i in 1:L
                out[i,j] *= n/(n+1)
                out[i,j] *= coef_tree[i,j,n+1]^2 / (n+1)
            end
        end
    end
end
