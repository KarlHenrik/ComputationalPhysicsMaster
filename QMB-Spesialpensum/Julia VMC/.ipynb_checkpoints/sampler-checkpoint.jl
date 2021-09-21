include("Wavefunctions/wavefunctions.jl")
include("hamiltonians.jl")

import Base.+, Base./

mutable struct Samples{T}
    E::Float64
    E2::Float64
    ∇Ψ::T
    ∇ΨE::T
end

Samples(wf::SimpleGaussian) = Samples(0.0, 0.0, 0.0, 0.0)
Samples(wf::Correlated)     = Samples(0.0, 0.0, 0.0, 0.0)

function sample!(samples, particles, wf::WaveFunction, ham::Hamiltonian)
    E = potential(particles, ham) + kinetic(particles, wf)
    ∇Ψ = paramDer(particles, wf)
    samples.E   += E
    samples.E2  += E^2
    samples.∇Ψ  += ∇Ψ
    samples.∇ΨE += ∇Ψ * E
    return
end

function gradient(samples)
    return 2 .* (samples.∇ΨE .- samples.E .* samples.∇Ψ)
end


function +(a::Samples, b::Samples)
    return Samples(a.E + b.E, a.E2 + b.E2, a.∇Ψ .+ b.∇Ψ, a.∇ΨE .+ b.∇ΨE)
end

function /(a::Samples, b::Number)
    return Samples(a.E/b, a.E2/b, a.∇Ψ/b, a.∇ΨE/b)
end
    
function /(b::Number, a::Samples)
    return Samples(a.E/b, a.E2/b, a.∇Ψ/b, a.∇ΨE/b)
end