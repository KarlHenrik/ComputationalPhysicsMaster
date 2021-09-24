import Random

abstract type Scheme end

include("statsschemes.jl")
include("optimizerschemes.jl")

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
        
    result = createResult(sampler_threads, dims, num) # GradientSampler/OneBodySampler: Weighted average
    return result                                     # BlockingSampler: Blocking method, then weighted average with error properly propegated
end