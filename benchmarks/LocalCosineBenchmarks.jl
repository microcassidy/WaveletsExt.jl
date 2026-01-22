module LocalCosineBenchmarks
using BenchmarkTools
import WaveletsExt:LocalCosineBasis
import WaveletsExt:LocalCosine

const SUITE = BenchmarkGroup(["localcosine"])
data = ones(1024)
analysis_g = addgroup!(SUITE, "analysis")
synthesis_g = addgroup!(SUITE, "sythesis")
bestbasis_g = addgroup!(SUITE, "bestbasis")


out = Array{Float64}(undef,1024,5,1)
lcb = LocalCosineBasis()
# cost = LocalCosineBasis.cost(results);
# f(d) = fit(d)

analysis_g["operator"] = @benchmarkable LocalCosine.analysis_operator!($out,$lcb,$data)
analysis_g["fold"] = @benchmarkable LocalCosine.fold!($lcb)

# bestbasis_g["cost"] =

end
