function plan_dct_iv(n,max_depth::Integer)::Vector{FFTWPlan}
    plans::Vector{FFTWPlan} = [FFTW.plan_r2r(zeros(n ÷ 2^j),FFTW.REDFT11) for j in 0:max_depth]
end

function coefficient_tree(n::Int,max_depth::Int)
    Matrix{Float64}(undef, max_depth + 1)
end

function reset!(tr::Matrix{Float64})
    fill!(tr,0)
end


function lcb_step!(h::Folded,out::Vector,x::Vector,F::FFTWPlan,bell::OrthonormalBell,j::Int,block_size::Int)
    # rpp = pseudopacket(blocks[end-block_size+1:end],bell,:right)
    left = pseudopacket(blocks[begin:block_size],bell,:left)
    @debug "lcb_step!"
    @views for i in 1:2^j-1
        centre = x[idx*block_size+begin:(idx+1)*block_size]
        right = x[idx*block_size+begin:(idx+1)*block_size]
        mul!(out[i*block_size+begin:(idx+1)*block_size],F,h)
        left = centre
    end
    centre = @view x[end-block_size+1:end]
    right = pseudopacket(centre,bell,:right)
    _fold!(h,centre,left,right,bell)
    mul!(view(out,length(out)-block_size+1:length(out)),F,h)
    nothing
end

function analysis(x::Vector{<:Real};max_depth::Int=4)
    L = length(x)
    J,n = dyad_length(x)
    m = div(n,1 << (max_depth + 1))
    bell = OrthonormalBell(m)

    plans = plan_dct_iv(n,max_depth)
    packet_space = prepare_packet_space(n,max_depth)
    coef_tree = Matrix{Float64}(undef, max_depth + 1)
    h = Folded{Float64}(L)
    function run(v::Vector = x)
        flush_packet_space!(packet_space)
        _lcb_step! = Base.Fix{4}((Base.Fix{3}(lcb_step!,v),bell))
        for j in 0:max_depth
            block_size = L ÷ 2^j
            _update_size!(h,block_size)
            _lcb_step!(h,view(coef_tree,:,j+1),plans[j+1],j,block_size)
        end
        return packet_space
    end
end
