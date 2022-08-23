mutable struct BlockingSampler <: Sampler
    savedEnergies::Vector{Float64}
    sampled_steps::Int64
    E_index::Int64
    BlockingSampler(sampled_steps) = new(Vector{Float64}(undef, sampled_steps), sampled_steps, 1)
end

function sample!(sampler::BlockingSampler, state, system)
    E = potential(state.positions, system.ham) + state.sample_m.kinetic
    sampler.savedEnergies[sampler.E_index] = E
    sampler.E_index += 1
    return sampler
end

struct Blocking <: StatsScheme end
createSampler(scheme::Blocking, wf, sampled_steps) = BlockingSampler(sampled_steps)



