module LocalCosine
using ..Utils
using FFTW
import FFTW:REDFT11,plan_r2r,r2rFFTWPlan
import LinearAlgebra:mul!,norm
const DCT_T = r2rFFTWPlan{Float64,Vector{Int32},false,1,Tuple{Int64}}

include("LocalCosine/utils.jl")
include("LocalCosine/cost.jl")
export LocalCosineBasis
include("LocalCosine/datastructure.jl")
include("LocalCosine/window.jl")
include("LocalCosine/fold.jl")
include("LocalCosine/unfold.jl")

export LocalCosineBasis, analysis_operator, analysis_operator!
include("LocalCosine/analysis.jl")
export LocalCosineBasis, synthesis_operator, synthesis_operator!
include("LocalCosine/synthesis.jl")
#foo
end
