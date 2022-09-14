struct OneBodySampler <: Sampler
    oneBodyDensity::Vector{Int64}
    start::Float64
    stop::Float64
    length::Int64
    step::Float64
    num::Int64
    sampled_steps::Int64
end
function OneBodySampler(start, stop, length, num, sampled_steps)
    return OneBodySampler(zeros(Int64, length), Float64(start), Float64(stop), length, (stop - start) / (length - 1), num, sampled_steps)
end

function sample!(sampler::OneBodySampler, walker, wf)
    (; positions) = walker
    for pos in positions
        distance = abs(pos)
        bin = 1
        while (distance > bin * sampler.step) && (bin <= sampler.length - 1)
            bin += 1
        end
        sampler.oneBodyDensity[bin] += 1
    end
    return sampler
end