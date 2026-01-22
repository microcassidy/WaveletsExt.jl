import DataStructures: Queue
Base.@kwdef mutable struct VTree{T<:Number} <: AbstractVector{T}
    data::AbstractVector{T}
    left::Union{VTree{T},Nothing} = nothing
    right::Union{VTree{T},Nothing} = nothing
    pruned::Bool = false
    VTree(val::AbstractArray{T}) where T <: Number = new{T}(val,nothing,nothing)
end
ssq(v::VTree{<:Number}) = sum(map(x-> abs(x)^2,v))
ssq(v::Nothing) = 0
Base.getindex(v::VTree{<:Number},i::Int) = v.data[i]
Base.length(v::VTree{<:Number},i::Int) = length(v.data)
Base.size(v::VTree{<:Number}) = size(v.data)
pruned(tr::VTree) = tr.pruned
pruned(t::Nothing) = true
function _setprune!(t::VTree,val::Bool)
    !isnothing(t.left) && setfield!(t.left,:pruned,val)
    !isnothing(t.right) && setfield!(t.right,:pruned,val)
    nothing
end



_gt(func::Function,tr::VTree{<:Number}) = func(tr) > func(tr.left) + func(tr.right)

isleaf(vt::VTree) = isnothing(vt.left) && isnothing(vt.right)

function prune(root::VTree{<:Number},pred::Function = Base.Fix1(_gt, ssq))
    """
    prune marker is inclusive of the tr and its children
    """
    isleaf(root) && return root
    queue = Queue{typeof(root)}()
    # prune_queue = Queue{typeof(root)}()
    push!(queue,root)
    while length(queue) > 0
        tr = popfirst!(queue)
        if pred(tr)
            _setprune!(tr,true)
        else
            push!(queue,tr.left)
            push!(queue,tr.right)
        end
    end
end



function _insert!(tr,val;left)::Bool
    !isnothing(ifelse(left,tr.left,tr.right)) && return false
    if left
        tr.left = VTree(val)
    else
        tr.right = VTree(val)
    end
    return true
end
function pushtr!(root::Union{Nothing,VTree{T}},
                   val::AbstractVector{T}) where T
    if isnothing(root)
        return VTree(val)
    end
    queue = Queue{Union{VTree{T},Nothing}}()
    push!(queue,root)
    while length(queue) > 0
        tr = popfirst!(queue)
        if _insert!(tr,val;left=true)
            return root
        else
            push!(queue,tr.left)
        end
        if _insert!(tr,val;left=false)
            return root
        else
            push!(queue,tr.right)
        end
    end
end

# import Printf
print_tr(vt::Nothing,prefix::String;is_last) = println(prefix, is_last ? "" : "│")
lrpr(tr) = (pruned(tr.left) || pruned(tr.right)) || isleaf(tr)
function print_tr(vt::VTree{<:Number},prefix::String;is_last::Bool,view=:pruned)
    if view == :pruned
        println(prefix, is_last ? "└──" : "├──", lrpr(vt) ? round(ssq(vt);digits=3) : "*┐")
    else
        println(prefix, is_last ? "└──" : "├──", round(ssq(vt);digits=3))
    end
end
has_right(vt::VTree) = !isnothing(vt.right)
print_tree(t::Nothing,prefix::String="",is_left::Bool=true;view) = nothing

function print_tree(tr::VTree, prefix::String="", is_left::Bool=true;view=:pruned)
    islast = isleaf(tr)

    IL = prefix * (!is_left ? "│   " : "    ")
    function _pruned()
        if pruned(tr)
            print_tr(nothing,prefix;is_last=is_left,view=view)
        else
            print_tr(tr,prefix;is_last=is_left,view=view)
            !pruned(tr.right) && print_tree(tr.right,IL,false;view=:pruned)
            !pruned(tr.left) && print_tree(tr.left,IL;view=:pruned)
        end
    end

    function full()
        print_tr(tr,prefix;is_last=is_left,view=view)
        print_tree(tr.right,IL,false;view=:full)
        print_tree(tr.left,IL;view=:full)
    end

    if view == :pruned
        _pruned()
    elseif view == :full
        full()
    else
        error("view $(view) not recognised")
    end
end

function window(η::Int, n::Int)
    β = ones(n)
    center = (n + 1) / 2
    @inline s(t) = sin(π * (t/1000) / (2π) )
    β[1:η] = map(t-> 1 + s(t), 0:η-1)
    β[end-(η-1):end] = reverse(β[1:η])
    return β
end


function modulate(M,j)
    L = exp2(-j) * M |> Int
    M2J = M*exp2(-j)
    ps = 0:2^j

    map(p-> exp(im * p * M2J),ps)
end

function dyadic_decomposition(height=10)
    MIN_BASE2 = 4
    M = exp2(height) |> Int
    plans = [plan_dct(zeros(1024 ÷ 2^j,2^j),1) for j in 0:MIN_BASE2]
    modulations = map(j -> modulate(M,j),0:MIN_BASE2)

    η = Int(exp2(MIN_BASE2-1))
    g = Base.Fix1(window,η)
    #I think that progressive dyadic decomposition already forms a tree. No neeed to organize
    function decompose(x::Vector{Float64})
        N = 2^height
        @assert length(x) == N
        tree = []
        ssqtree = Complex{Float64}[]
        root = nothing
        for j in 0:MIN_BASE2
            d = exp2(j) |> Integer
            L = N ÷ 2^j
            @debug "L:$L"
            ℱ = plans[j + 1]
            _window = repeat(g(L),1,d)
            partitions = reshape(x,(L ,2^j))
            @debug size(partitions),size(_window)
            transforms = Complex.(ℱ * partitions)
            _,nc = size(transforms)
            modulation = modulations[j+1]
            for i in 1:nc
                transforms[:,i] .*= modulation[i]
            end

            transforms .*= _window
            for col in eachcol(transforms)
                root = pushtr!(root,col)
            end
        end
        prune(root)
        return root
    end
    return decompose
end

function best_basis(tr::Nothing,index_set::Vector,basis_coeff::Vector,level,index) end
#TODO:maybe move vectors external
function best_basis(tr::VTree,index_set,basis_coeff,level=0,index=0)
    q = Queue{typeof(tr)}
    @inline _isleaf(tr) = isleaf(tr) || ( pruned(tr.left) && pruned(tr.right) )

    if _isleaf(tr)
        push!(basis_coeff,tr)
        push!(index_set,(level,index))
    else
        @debug index_set
        best_basis(tr.left,index_set,basis_coeff,level + 1,index)
        best_basis(tr.right,index_set,basis_coeff, level + 1,index + 1)
    end
end

function get_offsets(index_set)
    M= 2^10
    bottom = 0
    offsets = Int64[]
    for (j,p) in index_set
        offset = exp2(-j)*M*(p+1) |> Int
        push!(offsets,offset)
    end
    return offsets
end


using Plots
export test_time_freq
energy(v::AbstractVector{<:Number}) = real(v'v)
function test_time_freq()
    samples = load_single()[2][1,:]
    decompose = dyadic_decomposition()
    results = []
    for i in 0:1
        chunk = samples[1024*i + 1:1024(i+1)]
        tree = decompose(chunk)
        # print_tree(tree;view=:full)
        println("-----PRUNED-------")
        # print_tree(tree;view=:pruned)
        println("------END---------")
        index_set = []
        coefficients = []
        best_basis(tree,index_set,coefficients)
        push!(results,(; :index_set => index_set, :coeffients => coefficients))
        offsets = get_offsets(index_set)
        plt = plot(chunk)

        vline!(plt, offsets)
        savefig(plt, "local_cosine_chunk$i.png")

        plt2 = plot(dct(chunk))
        v = reduce(vcat, coefficients)
        @assert length(v) == 1024
        _dct = dct(chunk)
        pltdct = plot(abs.(_dct),title= "1024 point dct")
        pltlocal = plot(abs.(v),title = "local cosine basis")
        savefig(plot(pltdct,pltlocal), "local_cosine_compare_chunk$i.png")

    end
    return results
end
