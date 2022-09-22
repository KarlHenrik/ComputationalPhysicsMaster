function blocking(wf, ham, metro; nthreads=1)
    block_steps = metro.sample_steps
    d = log2(block_steps)
    if d%1 != 0
        d = Int(floor(d))
        #println("Length of data = $(block_steps) is not a power of 2, reducing to $(Int(round((2^d),digits=0))) steps")
        block_steps = 2^d
    end
    distr_steps = distribute_steps(block_steps, nthreads)
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
    std = Statistics.std(savedEnergies, corrected=false) / sqrt(length(savedEnergies))
    if std < 1e-16
        E_err = std
    else
        E_err = sqrt(block(savedEnergies))
    end
    
    return BlockingResult(E, E_err, std)
end

function createResult(sampler::BlockingSampler)
    E, E_err, std = block(sampler.savedEnergies)
    
    return BlockingResult(E, E_err, std)
end