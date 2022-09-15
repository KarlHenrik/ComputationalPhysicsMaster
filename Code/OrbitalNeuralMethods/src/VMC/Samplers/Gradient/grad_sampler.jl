mutable struct GradientSampler{T} <: Sampler
    # Storing computed values in case of the step being denied
    E::Float64
    E2::Float64
    ∇Ψ::T
    ∇ΨE::T
    # The running totals
    tot_E::Float64
    tot_E2::Float64
    tot_∇Ψ::T
    tot_∇ΨE::T
    const sampled_steps::Int64
end
GradientSampler(wf::Union{SimpleGaussian, Correlated}, sampled_steps) = GradientSampler(0.0, 0.0, 0.0, 0.0, sampled_steps)
#GradientSampler(wf::RBM, sampled_steps) = GradientSampler(0.0, 0.0, zero.((wf.a, wf.b, wf.W)), zero.((wf.a, wf.b, wf.W)), sampled_steps)

function update_sample!(sampler::GradientSampler, positions, wf, ham)
    sampler.E = kinetic(positions, wf) + potential(positions, ham)
    sampler.E2 = sampler.E^2
    
    sampler.∇Ψ = paramDer!(sampler.∇Ψ, positions, wf)
    sampler.∇ΨE = sampler.∇Ψ .* E
    
    return sampler
end

function sample!(sampler::GradientSampler{T}) where T
    (; E, E2, ∇Ψ, ∇ΨE) = sampler
    
    sampler.tot_E   += E
    sampler.tot_E2  += E2

    # TODO this will not work for the neural network. Need to make custom structs and operators
    sampler.tot_∇Ψ  = sampler.tot_∇Ψ .+ ∇Ψ
    sampler.tot_∇ΨE = sampler.tot_∇ΨE .+ ∇ΨE
    return sampler
end

struct GradientResult{T} <: Result
    E::Float64
    E_std::Float64
    E2::Float64
    gradient::T
end

function createResult(samplers::Vector{GradientSampler{T}}) where T
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    E = 0.0
    E2 = 0.0
    ∇Ψ = zero.(samplers[1].tot_∇Ψ)
    ∇ΨE = zero.(samplers[1].tot_∇ΨE)
    for sampler in samplers
        E +=         sampler.tot_E    / totalSteps
        E2 +=        sampler.tot_E2   / totalSteps
        ∇Ψ = ∇Ψ .+   sampler.tot_∇Ψ  ./ totalSteps
        ∇ΨE = ∇ΨE .+ sampler.tot_∇ΨE ./ totalSteps
    end
    
    E_std = sqrt(abs(E2 - E^2))
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end

function createResult(sampler::GradientSampler)
    E =   sampler.tot_E    / sampler.sampled_steps
    E2 =  sampler.tot_E2   / sampler.sampled_steps
    ∇Ψ =  sampler.tot_∇Ψ  ./ sampler.sampled_steps
    ∇ΨE = sampler.tot_∇ΨE ./ sampler.sampled_steps
    
    E_std = sqrt(E2 - E^2)
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end