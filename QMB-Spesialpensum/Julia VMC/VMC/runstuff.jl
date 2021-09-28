function run_scheme(wf, ham, metro, dimsnum, nthreads, scheme::StatsScheme)
    result = vmc(wf, ham, metro, dimsnum, nthreads, scheme)
    return result
end

function run_scheme(wf, ham, metro, dims, num, nthreads, scheme::OptimizerScheme)
    # Saving the sampled values and wavefunctions from each step to analyze later
    results = Vector{GradientResult}(undef, scheme.maxiter)
    wavefunctions = Vector{WaveFunction}(undef, scheme.maxiter)

    for i in 1:scheme.maxiter
        # Running the vmc calculation the get the gradient for this wavefunction
        result = vmc(wf, ham, metro, dims, num, nthreads, scheme)
        # Saving the results for this step for later analysis
        results[i] = result
        wavefunctions[i] = wf
        # Computing the gradient using a chosen optimization method
        grad = update_gradient(result.gradient, scheme) # This is where the optimization scheme is applied (Gradient Descent, Adam, etc.)
        # Stops the calculation if the gradient is smaller than the given tolerance
        if la.norm(grad) < scheme.tol
            results = results[1:i]
            wavefunctions = wavefunctions[1:i]
            break
        end
        # Applying gradient
        wf = applyGradient(wf, grad)

        # Output to keep track of progress
        if i % 5 == 0
            print("\rE = $(round(result.E, digits = 3)) iter = $i/$(scheme.maxiter)")
        end
    end
    
    return results, wavefunctions
end

function vmc(wf, ham, metro, dims, num, nthreads, scheme)
    steps_threads = [metro.sample_steps÷nthreads for i in 1:nthreads] # Dividing the steps equally between threads
    steps_threads[1] += metro.sample_steps % nthreads # The first thread also gets the steps that could not be evenly divided
    sampler_threads = [createSampler(scheme, wf, steps_threads[i]) for i in 1:nthreads] # Sampler for each thread
    
    Threads.@threads for i = 1:nthreads
        rng = Random.MersenneTwister()
        sampler = sampler_threads[i]
        particles = Particles(dims, num, rng, wf)
        
        for j in 1:metro.equil_steps
            metro_step!(particles, wf, rng, metro)
        end
        for j in 1:steps_threads[i]
            metro_step!(particles, wf, rng, metro)
            sample!(sampler, particles, wf, ham)
        end
    end
    
    result = createResult(sampler_threads) # GradientSampler/OneBodySampler: Weighted average
    return result                                     # BlockingSampler: Blocking method, then weighted average with error properly propegated
end