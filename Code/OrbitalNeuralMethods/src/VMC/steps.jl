function distribute_steps(steps, nthreads)
    distr_steps = [steps÷nthreads for i in 1:nthreads]
    distr_steps[1] += steps % nthreads
    return distr_steps
end

function steps!(samplers, wf_template, ham, metro)
    Threads.@threads for i = eachindex(samplers)
        # Move particles for a while to end up in a more likely state than the initial random state
        walker = Walker(wf_template)
        wf = private_wf(wf_template, walker.positions) # Creating a wf for this thread and these starting positions
        walker = equil_steps!(walker, wf, metro)

        sampler = samplers[i]
        update_sample!(sampler, walker, wf, ham)
        sampler, walker = sampled_steps!(sampler, walker, wf, ham, metro)
    end
    
    return samplers
end

function equil_steps!(walker, wf, metro)
    for i in 1:metro.equil_steps
        metro_step!(walker, wf, metro)
    end
    return walker
end

function sampled_steps!(sampler, walker, wf, ham, metro)
    for i in 1:sampler.sampled_steps
        accepted = metro_step!(walker, wf, metro)
        if accepted
            sampler = update_sample!(sampler, walker, wf, ham)
        end
        sampler = sample!(sampler)
    end
    return sampler, walker
end