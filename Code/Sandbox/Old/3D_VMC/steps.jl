function distribute_steps(steps, nthreads)
    distr_steps = [steps÷nthreads for i in 1:nthreads]
    distr_steps[1] += steps % nthreads
    return distr_steps
end

function steps!(samplers, wf, ham, metro)
    Threads.@threads for i = eachindex(samplers)
        # Move particles for a while to end up in a more likely state than the initial random state
        equil_walker = EquilWalker(wf, metro)
        wf_priv = private_wf(wf, equil_walker.positions) # Creating a wf for this thread and these starting positions
        equil_steps!(equil_walker, wf_priv, ham, metro)
        
        # Move particles and record the values of interest
        sampler = samplers[i]
        walker = SampledWalker(wf_priv, metro, sampler, equil_walker)
        sampler, walker = sampled_steps!(sampler, walker, wf_priv, ham, metro)
    end
    
    return samplers
end

function equil_steps!(walker, wf, ham, metro)
    for i in 1:metro.equil_steps
        walker = metro_step!(walker, wf, metro, ham)
    end
    return walker
end

function sampled_steps!(sampler, walker, wf, ham, metro)
    for i in 1:sampler.sampled_steps
        walker = metro_step!(walker, wf, metro, ham)
        sampler = sample!(sampler, walker, wf)
    end
    return sampler, walker
end