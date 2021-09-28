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

createSampler(scheme::OptimizerScheme, wf, sampled_steps) = GradientSampler(wf, sampled_steps)

struct GradientDescent <: OptimizerScheme
    lr::Float64
    maxiter::Int64
    tol::Float64
end
GradientDescent(;lr, maxiter, tol) = GradientDescent(lr, maxiter, tol)

function update_gradient(gradient, scheme::GradientDescent)
    return gradient .* scheme.lr
end

struct GradientResult{T} <: Result
    E::Float64
    E_std::Float64
    E2::Float64
    gradient::T
end

function createResult(sampler::GradientSampler)
    E = sampler.E / sampler.sampled_steps
    E2 = sampler.E2 / sampler.sampled_steps
    ∇Ψ = sampler.∇Ψ ./ sampler.sampled_steps
    ∇ΨE = sampler.∇ΨE ./ sampler.sampled_steps
    
    E_std = sqrt(E2 - E^2)
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end

function createResult(samplers::Vector{GradientSampler{T}}) where T
    E = 0.0
    E2 = 0.0
    ∇Ψ = zero(samplers[1].∇Ψ)
    ∇ΨE = zero(samplers[1].∇ΨE)
    
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    for sampler in samplers
        E += sampler.E / totalSteps
        E2 += sampler.E2 / totalSteps
        ∇Ψ += sampler.∇Ψ ./ totalSteps
        ∇ΨE += sampler.∇ΨE ./ totalSteps
    end
    E_std = sqrt(abs(E2 - E^2))
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end