function blocking(wf, ham, metro; nthreads=1)
    distr_steps = distribute_steps(metro.sample_steps, nthreads)
    samplers = [BlockingSampler(distr_steps[i]) for i in 1:nthreads]
    
    samplers = steps!(samplers, wf, ham, metro)
    
    block_result = createResult(samplers)
    return block_result
end

struct BlockingResult <: Result
    E::Float64
    E_err::Float64
    std::Float64
end

function createResult(samplers::Vector{BlockingSampler})
    savedEnergies = vcat([s.savedEnergies for s in samplers]...)
    
    E = Statistics.mean(savedEnergies)
    E_err = sqrt(block(savedEnergies))
    std = Statistics.std(savedEnergies, corrected=false) / sqrt(length(savedEnergies))
    
    return BlockingResult(E, E_err, std)
end

function createResult(sampler::BlockingSampler)
    E, E_err, std = block(sampler.savedEnergies)
    
    return BlockingResult(E, E_err, std)
end