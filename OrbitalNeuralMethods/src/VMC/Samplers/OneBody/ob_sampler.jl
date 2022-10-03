struct OneBodySampler <: Sampler
    oneBodyDensity_current::Vector{Int64}
    oneBodyDensity::Vector{Int64}
    start::Float64
    stop::Float64
    length::Int64
    step::Float64
    n::Int64
    sampled_steps::Int64
end
function OneBodySampler(start, stop, length, n, sampled_steps)
    return OneBodySampler(zeros(Int64, length), zeros(Int64, length),
                            Float64(start), Float64(stop), length, (stop - start) / (length - 1),
                            n, sampled_steps)
end

function update_sample!(sampler::OneBodySampler, walker::Walker, wf, ham)
    (; start, stop, step, length, oneBodyDensity_current) = sampler
    (; positions) = walker
    oneBodyDensity_current .= 0;
    @inbounds for pos in positions
        if pos < start
            oneBodyDensity_current[1] += 1
        elseif pos > stop
            oneBodyDensity_current[length] += 1
        else
            bin = 1
            while (pos > start + bin * step)
                bin += 1
            end
            oneBodyDensity_current[bin] += 1
        end
    end
    return sampler
end

function sample!(sampler::OneBodySampler)
    sampler.oneBodyDensity .+= sampler.oneBodyDensity_current
    return sampler
end