function blocking(wf, ham, metro, nthreads)
    distr_steps = distribute_steps(metro.sampled_steps, nthreads)
    samplers = [BlockingSampler(distr_steps[i]) for i in 1:nthreads]
    
    samplers = vmc!(samplers, wf, ham, metro)
    
    block_result = CreateResult(samplers)
    return block_result
end

struct BlockingResult <: Result
    E::Float64
    E_err::Float64
    std::Float64
end

function createResult(samplers::Vector{BlockingSampler})
    savedEnergies = vcat([s.savedEnergies for s in samplers]...)
    
    E, E_err, std = block(savedEnergies)
    
    return BlockingResult(E, E_err, std)
end

function createResult(sampler::BlockingSampler)
    E, E_err, std = block(sampler.savedEnergies)
    
    return BlockingResult(E, E_err, std)
end