module LocalCosine
import FFTW:REDFT11,plan_r2r,r2rFFTWPlan
using FFTW
# import FFTW:FFTWPlan
import LinearAlgebra:mul!
const DCT_T = r2rFFTWPlan{Float64,Vector{Int32},false,1,Tuple{Int64}}

include("utils.jl")
include("cost.jl")
export LocalCosineBasis
include("datastructure.jl")
include("window.jl")
include("fold.jl")

export LocalCosineBasis, analysis_operator, analysis_operator!
include("analysis.jl")
#foo
end
