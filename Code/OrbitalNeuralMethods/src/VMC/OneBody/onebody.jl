
function onebody(wf, ham, metro; start, stop, length, nthreads = 1)
    (; n) = wf
    distr_steps = distribute_steps(metro.sample_steps, nthreads)
    samplers = [OneBodySampler(start, stop, length, dims, n, distr_steps[i]) for i in 1:nthreads]
    
    samplers = steps!(samplers, wf, ham, metro)
    
    ob_result = createResult(samplers)
    return ob_result
end


struct OneBodyResult{T} <: Result
    oneBodyDensity::T
    radius::T
end

function createResult(sampler::OneBodySampler)
    (; oneBodyDensity, sample_steps, n, step, start, stop) = sampler
    return OneBodyResult(oneBodyDensity ./ sample_steps ./ n ./ step, [i for i in start:step:stop])
end

function createResult(samplers::Vector{OneBodySampler})
    (; length, n, step, start, stop) = samplers[1]
    oneBodyDensity = zeros(Float64, length)
    
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    for sampler in samplers
        oneBodyDensity .+= sampler.oneBodyDensity ./ totalSteps ./ n ./ step
    end
    
    return OneBodyResult(oneBodyDensity, [i for i in start:step:stop])
end
