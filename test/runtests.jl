ENV["GKSwstype"] = "100"
ENV["PLOTS_TEST"] = "true"

const TEST_GROUP = let g = get(ENV, "TEST_GROUP", nothing)
    g === nothing && length(ARGS) > 0 ? ARGS[1] : (g === nothing ? nothing : g)
end

using 
    Test,
    TestItems,
    TestItemRunner,
    Distributions,
    ImageQualityIndexes,
    Random,
    Plots,
    Statistics,
    Wavelets,
    WaveletsExt,
    SparseArrays


# @testitem  "Utils" begin include("utils.jl") end
# @testitem  "Transforms" begin include("transforms.jl") end
# @testitem  "Wavelet Multiplication" begin include("wavemult.jl") end
# @testitem  "Best Basis" begin include("bestbasis.jl") end
# @testitem  "Denoising" begin include("denoising.jl") end
# @testitem  "LDB" begin include("ldb.jl") end
# @testitem  "Visualizations" begin include("visualizations.jl") end

@testset  "LocalCosineBasis" begin include("localcosinebasis.jl") end
