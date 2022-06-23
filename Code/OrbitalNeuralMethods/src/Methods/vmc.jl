import LinearAlgebra as la
import StaticArrays as sa
import Random
import Statistics

# Independent
include("VMC/Wavefunctions/wavefunctions.jl") # Computes ratio, QF, parameter derivative, kinetic energy (WaveFunction)
include("VMC/hamiltonians.jl") # Computes potential energy (Hamiltonian)
# Needs WaveFunctions
include("VMC/particles.jl") # Saves particle positions (Particles)
include("VMC/metropolis.jl") # Moves the particles according to some metropolis algorithm (Metropolis/Importance)
# Needs WaveFunctions and Hamiltonians
include("VMC/Samplers/samplers.jl") # Chooses what is sampled during the calculation, and what is then computed after with the sampled values (GradientDescent/OneBody/Blocking <: Scheme)
                                # Samples values during the calculation (Gradient/OneBody/Blocking-Sampler)
                                # Computes results after (Gradient/OneBody/Blocking-Result)
# Needs everything
include("VMC/runstuff.jl")