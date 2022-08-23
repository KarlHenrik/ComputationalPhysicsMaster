function distribute_steps(steps, nthreads)
    distr_steps = [steps÷nthreads for i in 1:nthreads]
    distr_steps[1] += steps % nthreads
    return distr_steps
end

function steps!(samplers, wf, ham, metro)
    Threads.@threads for i = 1:length(samplers)
        # Move particles for a while to end up in a more likely state than the initial random state
        equil_walker = EquilWalker(wf, metro)
        equil_steps!(equil_walker, wf, ham, metro)
        starting_pos = equil_walker.positions
        
        # Move particles and record the values of interest
        sampler = samplers[i]
        walker = SampledWalker(starting_pos, wf, metro, sampler)
        sampler, walker = sampled_steps!(sampler, walker, wf, ham, metro)
    end
    
    return samplers
end

function equil_steps!(state, system, metro)
    for i in 1:metro.equil_steps
        walker = metro_step!(walker, wf, metro)
    end
    return walker
end

function sampled_steps!(sampler, walker, wf, ham, metro)
    for i in 1:sampler.sampled_steps
        walker = metro_step!(walker, wf, metro)
        sampler = sample!(sampler, walker, wf, ham)
    end
    return sampler, walker
end