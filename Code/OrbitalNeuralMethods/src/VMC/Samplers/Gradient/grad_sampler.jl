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
GradientSampler(wf::Union{SimpleGaussian, Correlated}, sampled_steps) = GradientSampler(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, sampled_steps)
#GradientSampler(wf::RBM, sampled_steps) = GradientSampler(0.0, 0.0, zero.((wf.a, wf.b, wf.W)), zero.((wf.a, wf.b, wf.W)), sampled_steps)
function GradientSampler(wf::SlaterNN, sampled_steps)
    ∇Ψ = paramDerHolder(wf.nn)
    ∇ΨE = paramDerHolder(wf.nn)
    tot_∇Ψ = paramDerHolder(wf.nn)
    tot_∇ΨE = paramDerHolder(wf.nn)
    GradientSampler(0.0, 0.0, ∇Ψ, ∇ΨE, 0.0, 0.0, tot_∇Ψ, tot_∇ΨE,sampled_steps)
end

function update_sample!(sampler::GradientSampler, walker, wf, ham)
    (; positions) = walker
    sampler.E = kinetic(positions, wf) + potential(positions, ham)
    sampler.E2 = sampler.E^2
    
    sampler.∇Ψ = paramDer!(sampler.∇Ψ, positions, wf)
    sampler.∇ΨE = setmul!(sampler.∇ΨE, sampler.∇Ψ, sampler.E)
    
    return sampler
end
setmul!(a, b, c) = b * c

function sample!(sampler::GradientSampler{T}) where T
    (; E, E2, ∇Ψ, ∇ΨE) = sampler
    
    sampler.tot_E   += E
    sampler.tot_E2  += E2
    sampler.tot_∇Ψ  = add!(sampler.tot_∇Ψ, ∇Ψ)
    sampler.tot_∇ΨE = add!(sampler.tot_∇ΨE, ∇ΨE)

    return sampler
end
add!(a, b) = a + b

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
    

    E = samplers[1].tot_E                      / totalSteps
    E2 = samplers[1].tot_E2                    / totalSteps
    ∇Ψ = samplers[1].tot_∇Ψ               ./ totalSteps
    ∇ΨE = samplers[1].tot_∇ΨE             ./ totalSteps

    
    for sampler in samplers[2:end]
        E +=         sampler.tot_E         / totalSteps
        E2 +=        sampler.tot_E2        / totalSteps

        setmul!.(sampler.tot_∇Ψ, sampler.tot_∇Ψ, 1 / totalSteps)
        ∇Ψ = add!.(∇Ψ, sampler.tot_∇Ψ)
        setmul!.(sampler.tot_∇ΨE, sampler.tot_∇ΨE, 1 / totalSteps)
        ∇ΨE = add!.(∇ΨE, sampler.tot_∇ΨE)
    end
    
    E_std = sqrt(abs(E2 - E^2))
    gradient = setmul!.(∇Ψ, ∇Ψ, -E)
    gradient = add!(gradient, ∇ΨE)
    gradient = setmul!.(gradient, gradient, 2)
    #gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
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