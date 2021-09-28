struct BlockingSampler <: Sampler
    savedEnergies::Vector{Float64}
    sampled_steps::Int64
    BlockingSampler(sampled_steps) = new(Vector{Float64}(undef, sampled_steps), sampled_steps)
end

function sample!(sampler::BlockingSampler, particles, wf::WaveFunction, ham::Hamiltonian)
    # TODO
    return
end

struct Blocking <: StatsScheme end
createSampler(scheme::Blocking, wf, sampled_steps) = BlockingSampler(sampled_steps)

struct BlockingResult <: Result
    E::Float64
    E_err::Float64
end

function createResult(sampler::BlockingSampler)
    return BlockingResult(E, E_err)
end

function createResult(samplers::Vector{BlockingSampler})
    s = samplers[1]
    E = 0.0
    E2 = 0.0
    gradient = zero(s.∇Ψ)
    
    for sampler in samplers
        E += sampler.E / sampler.sampled_steps
        E2 += sampler.E2 / sampler.sampled_steps
        gradient += gradient(sampler) / sampler.sampled_steps
    end
    
    return GradientResult(E, E2, gradient)
end