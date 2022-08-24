mutable struct GradientSampler{T} <: Sampler
    E::Float64
    E2::Float64
    ∇Ψ::T
    ∇ΨE::T
    const sampled_steps::Int64
end
GradientSampler(wf::SimpleGaussian, sampled_steps) = GradientSampler(0.0, 0.0, 0.0, 0.0, sampled_steps)
GradientSampler(wf::RBM, sampled_steps) = GradientSampler(0.0, 0.0, zero.((wf.a, wf.b, wf.W)), zero.((wf.a, wf.b, wf.W)), sampled_steps)

function sample!(sampler::GradientSampler, walker)
    (;positions, accepted, samp_muts) = walker
    
    E = samp_muts.E
    
    ∇Ψ = samp_muts.paramDer
    sampler.E   += E
    sampler.E2  += samp_muts.E2
    sampler.∇Ψ  = sampler.∇Ψ .+ ∇Ψ
    sampler.∇ΨE = sampler.∇ΨE .+ ∇Ψ .* E
    return sampler
end

createSampler(scheme::OptimizerScheme, wf, sampled_steps) = GradientSampler(wf, sampled_steps)


