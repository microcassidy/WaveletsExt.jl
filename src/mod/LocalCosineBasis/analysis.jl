function plan_dct_iv(n,max_depth::Integer)::Vector{FFTWPlan}
    plans::Vector{FFTWPlan} = [FFTW.plan_r2r(zeros(n ÷ 2^j),FFTW.REDFT11) for j in 0:max_depth]
end

function coefficient_tree(n::Int,max_depth::Int)
    Matrix{Float64}(undef, max_depth + 1)
end

function reset!(tr::Matrix{Float64})
    fill!(tr,0)
end


function lcb_step!(out::AbstractVector,
                   h::Folded,
                   F::FFTWPlan,
                   j::Int,
                   block_size::Int,
                   v::Vector,
                   bell::OrthonormalBell,
                   )
    endi = 2^j-1
    @views for idx in 0:endi
        left = idx == 0 ? pseudopacket(v[begin:block_size],bell,:left) : v[(idx-1)*block_size+begin:(idx)*block_size]
        centre = v[idx*block_size+begin:(idx+1)*block_size]
        right = idx == endi ? pseudopacket(centre,bell,:right) : v[idx*block_size+begin:(idx+1)*block_size]
        _fold!(h,centre,left,right,bell)
        out[idx*block_size+begin:(idx+1)*block_size] = F*h
    end
end

function cost(tr::Matrix{Float64})::Vector{Float64}
    nr,nc = size(tr)
    J = nc - 1
    ncost = (1 << nc) - 1
    cost_tr = zeros(ncost)
    i = 0
    for (j,c) in enumerate(eachcol(tr))
        nsplit = 2^(j-1)
        bsize = Int(nr / nsplit)
        @info "nsplit $nsplit bsize $bsize"
        @views for si in 0:nsplit-1
            v = c[si*bsize+1:(si+1)*bsize]
            cost_tr[i+begin]= v'v
            i+=1
        end
    end
    @assert i == ncost "ncost:$ncost i:$(i)"
    return cost_tr
end

flush_tree!(tr::AbstractArray{<:Number}) = fill!(tr,0)
function analysis(L::Int;max_depth::Int=4)
    J = Int(log2(L))
    m = div(L,1 << (max_depth + 1))
    bell = OrthonormalBell(m)

    plans = plan_dct_iv(L,max_depth)
    coef_tree = Matrix{Float64}(undef, (L, max_depth + 1))
    fill!(coef_tree,0)
    h = Folded(L)
    function run(v::Vector)
        flush_tree!(coef_tree)
        _lcb_step! = Base.Fix{6}(Base.Fix{6}(lcb_step!, v), bell)
        @assert length(v) == L
        flush_tree!(coef_tree)
        for j in 0:max_depth
            block_size = L ÷ 2^j
            _update_size!(h,block_size)
            _lcb_step!(view(coef_tree,:,j+1), h, plans[j+1], j, block_size)
        end
        return coef_tree
    end

    return run
end
