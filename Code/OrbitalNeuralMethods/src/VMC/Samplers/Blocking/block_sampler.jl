mutable struct BlockingSampler <: Sampler
    E::Float64
    E_index::Int64

    const savedEnergies::Vector{Float64}
    const sampled_steps::Int64

    function BlockingSampler(sampled_steps)
        return new(0.0, 1, Vector{Float64}(undef, sampled_steps), sampled_steps)
    end
end

function update_sample!(sampler::BlockingSampler, walker::Walker, wf, ham)
    sampler.E = kinetic(walker.positions, wf) + potential(walker.positions, ham)
    return sampler
end

function sample!(sampler::BlockingSampler)
    sampler.savedEnergies[sampler.E_index] = sampler.E
    sampler.E_index += 1
    return sampler
end


