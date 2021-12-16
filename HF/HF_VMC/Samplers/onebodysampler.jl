struct OneBodySampler{T} <: Sampler
    oneBodyDensity::T
    start::Float64
    stop::Float64
    length::Int64
    step::Float64
    dims::Int64
    num::Int64
    sampled_steps::Int64
end
OneBodySampler(start, stop, length, dims, num, sampled_steps) = OneBodySampler(zeros(Int64, length), start, stop, length, (stop - start) / (length - 1), dims, num, sampled_steps)

function sample!(sampler::OneBodySampler, particles, wf::WaveFunction, ham::Hamiltonian)
    for pos in particles.positions
        for dim in 1:sampler.dims
            distance = abs(pos[dim])
            bin = 1
            while (distance > bin * sampler.step) && (bin <= sampler.length - 1)
                bin += 1
            end
            sampler.oneBodyDensity[bin] += 1
        end
    end
    return
end

struct OneBody <: StatsScheme
    dims::Int64
    num::Int64
    start::Float64
    stop::Float64
    length::Int64
end
OneBody(dims, num; start, stop, length) = OneBody(dims, num, start, stop, length)
createSampler(scheme::OneBody, wf, sampled_steps) = OneBodySampler(scheme.start, scheme.stop, scheme.length, scheme.dims, scheme.num, sampled_steps)

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