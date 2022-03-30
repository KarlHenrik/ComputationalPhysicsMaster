function run_scheme(scheme::StatsScheme, system, metro, nthreads)
    result = vmc(scheme::StatsScheme, system, metro, nthreads)
    return result
end

function run_scheme(scheme::OptimizerScheme, system, metro, nthreads)
    # Saving the sampled values and wavefunctions from each step to analyze later
    results       = Vector{GradientResult}(undef, scheme.maxiter)
    wavefunctions = Vector{WaveFunction}(undef, scheme.maxiter)

    for i in 1:scheme.maxiter
        # Running the vmc calculation the get the gradient for this wavefunction
        result = vmc(scheme, system, metro, nthreads)
        # Saving the results for this step for later analysis
        results[i] = result
        wavefunctions[i] = system.wf
        # Computing the gradient using a chosen optimization method
        grad = update_gradient(result.gradient, scheme) # This is where the optimization scheme is applied (Gradient Descent, Adam, etc.)
        
        # Stops the calculation if the gradient is smaller than the given tolerance
        if la.norm(grad) < scheme.tol
            results = results[1:i]
            wavefunctions = wavefunctions[1:i]
            break
        end
        
        # Applying gradient
        new_wf = applyGradient(system.wf, grad)
        system = System(dims, num, new_wf, system.ham)

        # Output to keep track of progress
        print("\rE = $(round(result.E, digits = 3)) iter = $i/$(scheme.maxiter)                       ")
    end
    
    return results, wavefunctions
end

function vmc(scheme, system, metro, nthreads)
    steps_threads = [metro.sample_steps÷nthreads for i in 1:nthreads] # Dividing the steps equally between threads
    steps_threads[1] += metro.sample_steps % nthreads # The first thread also gets the steps that could not be evenly divided
    samplers = Vector{Sampler}(undef, nthreads)
    
    Threads.@threads for i = 1:nthreads
        state = EquilState(system, metro)
        # Move particles for a while to end up in a more likely state than the initial random state
        state = equil_steps!(state, system, metro)
        
        state = SampledState(state, system, metro, scheme)
        sampler = Sampler(scheme, system, steps_threads[i])
        # Move particles and record the values of interest
        state, sampler = sampled_steps!(state, sampler, system, metro)
        
        samplers[i] = sampler
    end
    
    result = createResult(samplers) # GradientSampler/OneBodySampler: Weighted average
    return result                   # BlockingSampler: Blocking method, then weighted average with error properly propegated
end

function equil_steps!(state, system, metro)
    for i in 1:metro.equil_steps
        state = metro_step!(state, system, metro)
    end
    return state
end

function sampled_steps!(state, sampler, system, metro)
    for i in 1:sampler.sampled_steps
        state = metro_step!(state, system, metro)
        sampler = sample!(sampler, state, system)
    end
    return state, sampler
end