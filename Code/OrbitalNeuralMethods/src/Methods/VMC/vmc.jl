# Sampler and Result setups
abstract type Sampler end
abstract type Result end

include("Blocking/block_algo.jl")
include("Blocking/block_sampler.jl")
include("Blocking/blocking.jl")

include("Gradient/grad_sampler.jl")
include("Gradient/optimize.jl")

include("OneBody/ob_sampler.jl")
include("OneBody/onebody.jl")

# 
include("Walker/muts.jl")
include("Walkers/walker.jl")
include("Walkers/consider.jl")

abstract type WaveFunction end
include("Wavefunctions/simplegaussian.jl")

include("hamiltonians.jl")
include("metropolis.jl")
include("steps.jl")
