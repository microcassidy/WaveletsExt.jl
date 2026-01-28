using WaveletsExt:LocalCosine
using WaveletsExt.LocalCosine:foldedge!,unfoldedge!,fold!,unfold!,LeftPacket,RightPacket,CentrePacket
import WaveletsExt.Utils:BinaryTree
const L = LocalCosine

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
# #@info ans
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
import WaveletsExt.BestBasis: bestbasis_treeselection


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



@testset "analysis/synthesis" begin
    # input = [sin(2π*t/1024) for t in 0:(1024*3)-1]
    input = ones(1024*3)
    output = similar(input)
    @assert length(input) == 3*1024
    # input = ones(Float64,1024*3)
    lcb = LocalCosineBasis()
    coef_tree = Array{Float64,3}(undef,1024,5,3)
    #@info "input: $(size(coef_tree))"
    analysis_operator!(coef_tree,lcb,input)
    #@info coef_tree[1:10,1:4,1]

    _cost = Vector{Float64}(undef, 2^(lcb.max_depth+1)- 1)
    #@info "cost len $(length(_cost))"
    fill!(_cost,0)
    # (i -> #@info coef_tree[1:10,i]).(1:3)
    cost!(_cost,coef_tree)
    #@info _cost
    # #@info _cost

    # # bestbasis_treeselection(costs::AbstractVector{T},
    # bb = bestbasis_treeselection(_cost,1024, :max)
    bb = BitVector(zeros(2^(lcb.max_depth+1)-1))
    fill!(bb,false)
    # bb[2:3] .= 1
    bb[1] = true
    # @assert sum(bb) == 2
    tm = bestbasis_treemask(BinaryTree, bb,1024,coef_tree)
    bb_decomp = coef_tree[tm]
    plan = SynthesisPlan(bb,1024)

    v = coef_tree[tm]
    cr = 1:5
    #@info tm[:,2,:]
    # @test all(tm[:,2,:] .== true)
    # @test all(tm[:,[1;3:size(tm,2)],:] .== false)
    synthesis_operator!(output, v, lcb,plan,bb)
    err = log10.((output - input).^2)
    plt = plot(err,label="err")
    vline!(plt,(1:6).*512,col=:red)

    savefig(plt,"test_synthesis.pdf")
    @test isapprox(output-input,0;atol=10^-5)
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
