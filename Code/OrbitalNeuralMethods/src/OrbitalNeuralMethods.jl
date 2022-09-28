module OrbitalNeuralMethods

# ------------ Exports ---------------

# Orbital systems
export System, Basis, reference_energy, particle_density
export SpatialSystem, SpatialBasis, HOBasis, SpinBasis
export PairingSystem, pairing_exact, pairing_MBPT2, Pairing
export CalogeroSutherland, HOCoulomb, HarmonicOscillator

# Orbital system solvers
export energy, corr_energy, update!, compute_ground_state!
export HFState, HF
export RHFState, RHF
export CCDState, CCD
export CCSDState, CCSD
export Alpha, DIIS

# VMC
export SimpleGaussian, Correlated, Slater, SlaterNN
export NeuralNetwork, Dense, Exp, Tanh, Sigmoid
export Metropolis, Importance
export blocking, onebody, optimize, GradientDescent, ADAM

# ------------ Imports ---------------

import LinearAlgebra as la
import Random
import Statistics
using TensorOperations: @tensor
import ForwardDiff as fd
import ReverseDiff as rd

# ------------ Source Files ---------------

include("Systems/system.jl")

include("Methods/mixer.jl")
include("Methods/util.jl")

include("Methods/hf.jl")
include("Methods/rhf.jl")
include("Methods/ccd.jl")
include("Methods/ccsd.jl")

include("VMC/vmc.jl")

end
