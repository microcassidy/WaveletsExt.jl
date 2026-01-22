using WaveletsExt:LocalCosine

data = ones(1024)
lcb = LocalCosineBasis()
results = analysis_operator(lcb,data)

fill!(lcb.packet.left,1)
fill!(lcb.packet.right,1)
fill!(lcb.packet.centre,1)
LocalCosine.fold!(lcb)
using Plots
savefig(plot(lcb.packet.centre), "unit_response_fold.png")

# using BenchmarkTools
# f() = LocalCosine.fold!(lcb)
# f()
# @benchmark f()
