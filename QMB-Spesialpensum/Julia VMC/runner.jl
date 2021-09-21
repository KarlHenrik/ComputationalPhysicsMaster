include("particles.jl")
include("Wavefunctions/wavefunctions.jl")
include("sampler.jl")
include("metropolis.jl")

import Random

function gradientDescent(wf, ham, metro, dimsnum, nthreads, lr, maxiter, tol)
    for i in 1:maxiter
        samples = run(wf, ham, metro, dimsnum, nthreads)
        grad = gradient(samples)
    
        wf = applyGradient(wf, grad .* lr)
        if i % 1 == 0
            print("\ralpha = $(round(wf.alpha, digits = 3)) iter = $i/$maxiter")
        end
        if abs(grad) < tol
            return wf
        end
    end
    return wf
end


function run(wf, ham, metro, dimsnum, nthreads)
    rng_threads = [Random.MersenneTwister() for i in 1:nthreads]
    samples_threads = [Samples(wf) for i in 1:nthreads]
    steps_threads = [metro.sample_steps÷nthreads for i in 1:nthreads]
    steps_threads[1] += metro.sample_steps % nthreads

    Threads.@threads for i = 1:nthreads
        rng = rng_threads[i]
        samples = samples_threads[i]
        particles = Particles(dimsnum..., rng, wf)
        
        for j in 1:metro.equil_steps
            metro_step!(particles, wf, rng, metro)
        end
        for j in 1:steps_threads[i]
            metro_step!(particles, wf, rng, metro)
            sample!(samples, particles, wf, ham)
        end
    end
    
    total_samples = Samples(wf)
    for samples in samples_threads
        total_samples += samples
    end
    return total_samples / metro.sample_steps
end;