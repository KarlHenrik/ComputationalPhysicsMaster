mutable struct BlockingSampler <: Sampler
    const savedEnergies::Vector{Float64}
    const sampled_steps::Int64
    E_index::Int64
    function BlockingSampler(sampled_steps)
        return new(Vector{Float64}(undef, sampled_steps), sampled_steps, 1)
    end
end

function sample!(sampler::BlockingSampler, walker, wf)
    E = walker.samp_muts.E
    sampler.savedEnergies[sampler.E_index] = E
    sampler.E_index += 1
    return sampler
end


