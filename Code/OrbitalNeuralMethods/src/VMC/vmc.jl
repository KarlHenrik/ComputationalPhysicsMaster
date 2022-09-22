include("Wavefunctions/Tools/fastDeterminant.jl")
include("Wavefunctions/Tools/layer.jl")
include("Wavefunctions/Tools/minisort.jl")

abstract type WaveFunction end
include("Wavefunctions/Gaussians/gaussian_common.jl")
include("Wavefunctions/Gaussians/simplegaussian.jl")
include("Wavefunctions/Gaussians/correlated.jl")

include("Wavefunctions/slater.jl")
include("Wavefunctions/neuralnetwork.jl")
include("Wavefunctions/slaterNN.jl")

include("walker.jl")
include("hamiltonians.jl")
include("metropolis.jl")
include("steps.jl")

# Sampler and Result setups
abstract type Sampler end
abstract type Result end

include("Samplers/Blocking/block_algo.jl")
include("Samplers/Blocking/block_sampler.jl")
include("Samplers/Blocking/blocking.jl")

include("Samplers/Gradient/grad_sampler.jl")
include("Samplers/Gradient/optimize.jl")

include("Samplers/OneBody/ob_sampler.jl")
include("Samplers/OneBody/onebody.jl")




