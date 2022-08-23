function onebody(wf, ham, metro; start, stop, length)
    dims, num = wf.dims, wf.num
    distr_steps = distribute_steps(metro.sampled_steps, nthreads)
    samplers = [OneBodySampler(start, stop, length, dims, num, distr_steps[i]) for i in 1:nthreads]
    
    samplers = vmc!(samplers, wf, ham, metro)
    
    ob_result = CreateResult(samplers)
    return ob_result
end

struct OneBodyResult{T} <: Result
    oneBodyDensity::T
    radius::T
end

function createResult(sampler::OneBodySampler)
    s = sampler
    return OneBodyResult(s.oneBodyDensity ./ s.sampled_steps ./ s.num ./ s.dims ./ s.step, [i for i in s.start:s.step:s.stop])
end

function createResult(samplers::Vector{OneBodySampler{T}}) where T
    s = samplers[1]
    oneBodyDensity = zeros(Float64, s.length)
    
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    for sampler in samplers
        oneBodyDensity .+= sampler.oneBodyDensity ./ totalSteps ./ s.num ./ s.dims ./ s.step
    end
    
    return OneBodyResult(oneBodyDensity, [i for i in s.start:s.step:s.stop])
end