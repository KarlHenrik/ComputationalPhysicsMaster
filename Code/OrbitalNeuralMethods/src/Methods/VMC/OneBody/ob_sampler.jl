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

function sample!(sampler::OneBodySampler, state, system)
    for pos in state.positions
        for dim in 1:sampler.dims
            distance = abs(pos[dim])
            bin = 1
            while (distance > bin * sampler.step) && (bin <= sampler.length - 1)
                bin += 1
            end
            sampler.oneBodyDensity[bin] += 1
        end
    end
    return sampler
end