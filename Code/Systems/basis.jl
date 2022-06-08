abstract type Basis end

abstract type SpatialBasis <: Basis end

include("BasisFunctions/HObasis.jl")
include("BasisFunctions/spinbasis.jl")
include("BasisFunctions/pairingbasis.jl")
