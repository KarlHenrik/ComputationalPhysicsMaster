include("wavefunctions.jl")
include("hamiltonians.jl")

import Base.+

mutable struct Samples{T<:Union{Float64, Array{Float64}}}
    E::Float64
    E2::Float64
    ∇Ψ::T # this will in general have any shape
    ∇ΨE::T # this will in general have any shape
end

Samples() = Samples(0.0, 0.0, 0.0, 0.0)

function sample!(samples, positions, wf::WaveFunction, ham::Hamiltonian, temp)
    E = potential(positions, ham, temp) + kinetic(positions, wf, temp)
    ∇Ψ = paramDer(positions, wf, temp)
    samples.E   += E
    samples.E2  += E^2
    samples.∇Ψ  += ∇Ψ
    samples.∇ΨE += ∇Ψ * E
    return
end


function +(a::Samples, b::Samples)
    return Samples(a.E + b.E, a.E2 + b.E2, a.∇Ψ .+ b.∇Ψ, a.∇ΨE .+ b.∇ΨE)
end