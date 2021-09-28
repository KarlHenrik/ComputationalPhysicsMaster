import LinearAlgebra as la
import StaticArrays as sa
import Random

# Independent
include("Wavefunctions/wavefunctions.jl") # Computes ratio, QF, parameter derivative, kinetic energy (WaveFunction)
include("hamiltonians.jl") # Computes potential energy (Hamiltonian)
# Needs WaveFunctions
include("particles.jl") # Saves particle positions (Particles)
include("metropolis.jl") # Moves the particles according to some metropolis algorithm (Metropolis/Importance)
# Needs WaveFunctions and Hamiltonians
include("Samplers/samplers.jl") # Chooses what is sampled during the calculation, and what is then computed after with the sampled values (GradientDescent/OneBody/Blocking <: Scheme)
                                # Samples values during the calculation (Gradient/OneBody/Blocking-Sampler)
                                # Computes results after (Gradient/OneBody/Blocking-Result)
# Needs everything
include("runstuff.jl")