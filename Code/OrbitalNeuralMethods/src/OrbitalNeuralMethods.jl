module OrbitalNeuralMethods

# ------------ Exports ---------------

# Orbital systems
export System, Basis, reference_energy
export SpatialSystem, SpatialBasis, HOBasis, SpinBasis
export PairingSystem, pairing_exact, pairing_MBPT2, Pairing
export Interaction, CalogeroSutherland, ShieldedCoulomb, Coulomb, NonInteracting

# Orbital system solvers
export energy, corr_energy, update!, compute_ground_state!
export HFState, setup_HF
export RHFState, setup_RHF
export CCDState, setup_CCD
export CCSDState, setup_CCSD
export Alpha, DIIS

# VMC
export HarmonicOscillator
export SimpleGaussian, Correlated
export Metropolis, Importance
export blocking, onebody, optimize, GradientDescent

# ------------ Imports ---------------

import LinearAlgebra as la
import StaticArrays as sa
import Random
import Statistics
using TensorOperations: @tensor
import ForwardDiff as fd
import ReverseDiff as rd
import LinearAlgebra as la

# ------------ Source Files ---------------

include("Systems/system.jl")

include("Methods/mixer.jl")
include("Methods/util.jl")

include("Methods/hf.jl")
include("Methods/rhf.jl")
include("Methods/ccd.jl")
include("Methods/ccsd.jl")

include("Methods/VMC/vmc.jl")

end
