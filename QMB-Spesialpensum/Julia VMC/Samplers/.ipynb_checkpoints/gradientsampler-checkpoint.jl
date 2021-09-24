mutable struct GradientSampler{T} <: Sampler
    E::Float64
    E2::Float64
    ∇Ψ::T
    ∇ΨE::T
    sampled_steps::Int64
end
GradientSampler(wf::SimpleGaussian, sampled_steps) = GradientSampler(0.0, 0.0, 0.0, 0.0, sampled_steps)
GradientSampler(wf::Correlated, sampled_steps)     = GradientSampler(0.0, 0.0, 0.0, 0.0, sampled_steps)

function sample!(sampler::GradientSampler, particles, wf::WaveFunction, ham::Hamiltonian)
    E = potential(particles, ham) + kinetic(particles, wf)
    ∇Ψ = paramDer(particles, wf)
    sampler.E   += E
    sampler.E2  += E^2
    sampler.∇Ψ  += ∇Ψ
    sampler.∇ΨE += ∇Ψ * E
    return
end

function computeGradient(sampler::GradientSampler)
    return 2 .* (sampler.∇ΨE .- sampler.E .* sampler.∇Ψ)
end

struct GradientResult{T} <: Result
    E::Float64
    E_std::Float64
    E2::Float64
    gradient::T
end

function createResult(sampler::GradientSampler, dims, num)
    norm_fac = 1 / num / length(samplers)
    E = sampler.E / sampler.sampled_steps * norm_fac
    E2 = sampler.E2 / sampler.sampled_steps * norm_fac
    gradient = computeGradient(sampler) ./ sampler.sampled_steps .* norm_fac
    
    return GradientResult(E, sqrt(E2 - E^2), E2, gradient)
end

function createResult(samplers::Vector{GradientSampler{T}}, dims, num) where T
    s = samplers[1]
    E = 0.0
    E2 = 0.0
    gradient = zero(s.∇Ψ)
    
    norm_fac = 1 / num / length(samplers)
    for sampler in samplers
        E += sampler.E / sampler.sampled_steps * norm_fac
        E2 += sampler.E2 / sampler.sampled_steps * norm_fac
        gradient += computeGradient(sampler) ./ sampler.sampled_steps .* norm_fac
    end
    
    
    return GradientResult(E, sqrt(E2 - E^2), E2, gradient)
end