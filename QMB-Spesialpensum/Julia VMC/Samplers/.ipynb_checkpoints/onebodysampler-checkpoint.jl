mutable struct OneBodySampler{T} <: Sampler
    oneBodyDensity::T
    start::Int64
    stop::Int64
    length::Int64
    sampled_steps::Int64
end
OneBodySampler(start, stop, length, sampled_steps) = OneBodySampler([i for i in range(start, stop, length = length)], start, stop, length, sampled_steps)

function sample!(sampler::OneBodySampler, particles, wf::WaveFunction, ham::Hamiltonian)
    # TODO
    return
end

struct OneBodyResult{T} <: Result
    oneBodyDensity::T
end

function createResult(sampler::OneBodySampler, dims, num)
    
    return OneBodyResult()
end

function createResult(samplers::Vector{OneBodySampler}, dims, num)
    oneBodyDensity = zero(samplers[1].oneBodyDensity)
    
    for sampler in samplers
        oneBodyDensity .+= sampler.oneBodyDensity ./ sampler.sampled_steps
    end
    
    return OneBodyResult(oneBodyDensity)
end