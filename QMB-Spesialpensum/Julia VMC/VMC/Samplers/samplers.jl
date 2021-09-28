abstract type Scheme end
abstract type StatsScheme <: Scheme end
abstract type OptimizerScheme <: Scheme end

abstract type Sampler end
abstract type Result end

include("blockingsampler.jl")
include("gradientsampler.jl")
include("onebodysampler.jl")