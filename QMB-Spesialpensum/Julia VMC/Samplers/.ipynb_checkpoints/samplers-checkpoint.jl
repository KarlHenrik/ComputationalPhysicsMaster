abstract type Sampler end
abstract type Result end

include("blockingsampler.jl")
include("gradientsampler.jl")
include("onebodysampler.jl")