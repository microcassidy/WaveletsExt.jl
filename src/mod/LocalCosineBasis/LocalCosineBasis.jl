module LocalCosineBasis
using FFTW
import FFTW:FFTWPlan
import LinearAlgebra:mul!

include("utils.jl")
include("datastructure.jl")
include("window.jl")
include("fold.jl")
include("analysis.jl")
end
