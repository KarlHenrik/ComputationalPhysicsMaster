#=
function onebody(wf, ham, metro; start, stop, length, nthreads = 1)
    (; dims, num) = wf
    distr_steps = distribute_steps(metro.sample_steps, nthreads)
    samplers = [OneBodySampler(start, stop, length, dims, num, distr_steps[i]) for i in 1:nthreads]
    
    samplers = steps!(samplers, wf, ham, metro)
    
    ob_result = CreateResult(samplers)
    return ob_result
end
=#
#=
struct OneBodyResult{T} <: Result
    oneBodyDensity::T
    radius::T
end

function createResult(sampler::OneBodySampler)
    (; oneBodyDensity, sample_steps, num, dims, step, start, stop) = sampler
    return OneBodyResult(oneBodyDensity ./ sample_steps ./ num ./ dims ./ step, [i for i in start:step:stop])
end

function createResult(samplers::Vector{OneBodySampler})
    (; length, num, dims, step) = samplers[1]
    oneBodyDensity = zeros(Float64, length)
    
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    for sampler in samplers
        oneBodyDensity .+= sampler.oneBodyDensity ./ totalSteps ./ num ./ dims ./ step
    end
    
    return OneBodyResult(oneBodyDensity, [i for i in s.start:s.step:s.stop])
end
=#