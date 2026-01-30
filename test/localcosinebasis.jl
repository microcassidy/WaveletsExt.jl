using WaveletsExt:LocalCosine
using WaveletsExt.LocalCosine:foldedge!,unfoldedge!,fold!,unfold!,LeftPacket,RightPacket,CentrePacket
import WaveletsExt.Utils:BinaryTree
import WaveletsExt.BestBasis: bestbasis_treeselection
const L = LocalCosine
using WaveletsExt:Visualizations




# import Test:

# data = ones(1024)
# lcb = LocalCosineBasis()
# # results = analysis_operator(lcb,data)

# fill!(lcb.packet.left,1)
# fill!(lcb.packet.right,1)
# fill!(lcb.packet.centre,1)
# LocalCosine.fold!(lcb)
# using Plots
# savefig(plot(lcb.packet.centre), "unit_response_fold.png")

# # using BenchmarkTools
# # f() = LocalCosine.fold!(lcb)
# # f()
# # @benchmark f()
# output = zeros(Float64,31)
# input = ones(Float64,1024,5,6)
# ans = reduce(vcat,[repeat([div(1024,2^j)],2^j) for j in 0:4])
# LocalCosine.cost!(output,input)
# @test ans ≈ output


@testset "folding" begin
    folds(l::LocalCosineBasis)=([l.packet.left ; zeros(lcb.N) ; l.packet.right])
    input = rand(Float64, 1024)

    lcb = LocalCosineBasis()
    fill!(lcb.packet.left,0)
    fill!(lcb.packet.right,0)
    copyto!(lcb.packet.centre,input)

    foldedge!(LeftPacket,lcb)
    foldedge!(RightPacket,lcb)

    unfold!(lcb)
    unfoldedge!(LeftPacket,lcb)
    unfoldedge!(RightPacket,lcb)
    @test lcb.packet.centre ≈ input
end
# using WaveletsExtBestBasis


# column_index(bit_index::Int) = log2(bit_index) |> floor |> Float64 |> x->x+1
# """
# convert an index to a row-range and a column
# """
# function row_range(bit_index::Int,signal_length::int)
#     @assert 2^floor(log2(signal_length)) == 2^log2(signal_legnth) "signal length must be dyadic"
#     #the column index is equivalent to the jth level decomposition
#     # the rows of level j are broken up into 2^j portions we just need to figure out where it is
#     idx == 1 && return 1:signal_length
#     parent_range = row_range(getparentindex(bit_index, :binary))
# end
# function rand_bestbasis!(bb,N,idx,level=0)
#    split = rand(Bool)
#    level > 4
#    v
#    if !split || level == 4
#        fill!(,true)
#        return
#    else
#        mp = div(rng.start + rng.stop,2)
#        split(bb[begin:2^level])
#        R = mp+1:
#     end

# end
@testset "LocalCosineBasis" begin
    @test_throws ArgumentError LocalCosineBasis(1023)
end
@testset "analysis/synthesis" begin
    function setup(l,callback::Function=(_->nothing))
        input = rand(l)
        output = similar(input)
        lcb = LocalCosineBasis()
        coef_tree = Array{Float64,3}(undef,1024,5,div(l,1024))
        _cost = zeros(Float64,2^(lcb.max_depth+1)- 1)
        return (; :tr=> coef_tree,:lcb=>lcb,:in=>input,:out=>output,:cost=>_cost)
    end
    BB(d) = bestbasis_treeselection(d.cost,1024, :max)
    cost(d) = cost!(d.cost,d.tr)
    Ψ(d::NamedTuple) = analysis_operator!(d.tr,d.lcb,d.in)

    Ψ⁻¹(d::NamedTuple,v::Vector{Float64},bb::BitVector) = synthesis_operator!(d.out,v,d.lcb,bb)

    even = setup(1024*4)


    let
        odd = setup(1024*3)
        Ψ(odd)
        cost(odd)
        bb = trues(length(odd.cost))
        prune!(bb,odd.cost,:max)

        # @info odd.cost[1:7]
        # return
        # bb = BB(odd)
        @info bb
        tm = bestbasis_treemask(BinaryTree,bb,1024,odd.tr)
        v = odd.tr[tm]
        Ψ⁻¹(odd,v,bb)
        @test odd.out ≈ odd.in
    end

    let
        even = setup(1024*4)
        Ψ(even)
        cost(even)
        bb = trues(length(even.cost))
        prune!(bb,even.cost,:max)

        # @info even.cost[1:7]
        # return
        # bb = BB(even)
        @info bb
        tm = bestbasis_treemask(BinaryTree,bb,1024,even.tr)
        v = even.tr[tm]
        Ψ⁻¹(even,v,bb)
        @test even.out ≈ even.in
    end


    let
        include("ecg_sig.jl") #returns a var data
        ecg = setup(1024)
        copyto!(ecg.in,ecg_data)
        Ψ(ecg)
        cost(ecg)
        bb = trues(length(ecg.cost))
        prune!(bb,ecg.cost,:max)
        p1 = plot_tfbdry(bb)
        savefig(p1,"best_basis_ecg.pdf")
        @info bb
        tm = bestbasis_treemask(BinaryTree,bb,1024,ecg.tr)
        v = ecg.tr[tm]
        Ψ⁻¹(ecg,v,bb)
        @test ecg.out ≈ ecg.in
    end






    # cost!(_cost,coef_tree)
    # bb = bestbasis_treeselection(_cost,1024, :max)
    # tm = bestbasis_treemask(BinaryTree, bb,1024,coef_tree)
    # bb_decomp = coef_tree[tm]
    # v = coef_tree[tm]
    # cr = 1:5
    # synthesis_operator!(output, v, lcb,bb)
    # err = log10.((output - input).^2)
    # plt = plot(err,label="err")
    # vline!(plt,(1:6).*512,col=:red)
    # savefig(plt,"test_synthesis.pdf")
    # @test input ≈ output
end

@testset "rowindexes" begin
    xs = 1:(2^5-1)
    for x in xs
        rng = getrowrange(BinaryTree,1024,x)
        @test iseven(length(rng))
    end
end

    # function SynthesisPlan(tr::BitVector,N::Int64)
# function analysis_operator!(coef_tree::Array{Float64,3},
