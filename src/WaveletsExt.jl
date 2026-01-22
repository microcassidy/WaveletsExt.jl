__precompile__()

module WaveletsExt

include("mod/Utils.jl")
include("mod/DWT.jl")
include("mod/ACWT.jl")
include("mod/BestBasis.jl")
include("mod/SIWT.jl")
include("mod/SWT.jl")
include("mod/Denoising.jl")
include("mod/LDB.jl")
include("mod/Visualizations.jl")
include("mod/WaveMult.jl")
include("mod/LocalCosine/LocalCosine.jl")

using Reexport
@reexport using .DWT,
                .BestBasis,
                .Denoising,
                .Utils,
                .LDB,
                .SWT,
                .SIWT,
                .ACWT,
                .Visualizations,
                .WaveMult,
                .LocalCosine

using BenchmarkTools

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.15
BenchmarkTools.DEFAULT_PARAMETERS.memory_tolerance = 0.01

const PARAMS_PATH = joinpath(dirname(@__FILE__), "..", "etc", "params.json")
const SUITE = BenchmarkGroup()
const MODULES = Dict("LocalCosine" => :LocalCosineBenchmarks)

load!(id::AbstractString; kwargs...) = load!(SUITE, id; kwargs...)

function load!(group::BenchmarkGroup, id::AbstractString; tune::Bool = false)
    modsym = MODULES[id]
    # modpath = joinpath(dirname(@__FILE__), id, "../benchmarks/$(modsym).jl")
    modpath = joinpath(dirname(@__FILE__), "..", "benchmarks/$(modsym).jl")
    Core.eval(WaveletsExt, :(include($modpath)))
    mod = Core.eval(WaveletsExt, modsym)
    modsuite = @invokelatest getglobal(mod, :SUITE)
    group[id] = modsuite
    if tune
        results = BenchmarkTools.load(PARAMS_PATH)[1]
        haskey(results, id) && loadparams!(modsuite, results[id], :evals)
    end
    return group
end

loadall!(; kwargs...) = loadall!(SUITE; kwargs...)

function loadall!(group::BenchmarkGroup; verbose::Bool = true, tune::Bool = false)
    for id in keys(MODULES)
        if verbose
            print("loading group $(repr(id))... ")
            time = @elapsed load!(group, id, tune=false)
            println("done (took $time seconds)")
        else
            load!(group, id, tune=false)
        end
    end
    if tune
        results = BenchmarkTools.load(PARAMS_PATH)[1]
        for (id, suite) in group
            haskey(results, id) && loadparams!(suite, results[id], :evals)
        end
    end
    return group
end

end
