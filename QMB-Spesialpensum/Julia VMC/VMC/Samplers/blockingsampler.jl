mutable struct BlockingSampler <: Sampler
    savedEnergies::Vector{Float64}
    sampled_steps::Int64
    E_index::Int64
    BlockingSampler(sampled_steps) = new(Vector{Float64}(undef, sampled_steps), sampled_steps, 1)
end

function sample!(sampler::BlockingSampler, particles, wf::WaveFunction, ham::Hamiltonian)
    E = potential(particles, ham) + kinetic(particles, wf)
    sampler.savedEnergies[sampler.E_index] = E
    sampler.E_index += 1
    return
end

struct Blocking <: StatsScheme end
createSampler(scheme::Blocking, wf, sampled_steps) = BlockingSampler(sampled_steps)

struct BlockingResult <: Result
    E::Float64
    E_err::Float64
    std::Float64
end

function createResult(sampler::BlockingSampler)
    E, E_err, std = block(sampler.savedEnergies)
    
    return BlockingResult(E, E_err, std)
end

function createResult(samplers::Vector{BlockingSampler})
    savedEnergies = vcat([s.savedEnergies for s in samplers]...)
    
    E, E_err, std = block(savedEnergies)
    
    return BlockingResult(E, E_err, std)
end

function block(x)
    n0 = length(x)
    n = length(x)
    d = log2(n)
    d_int = Int(floor(d))
    s = zeros(d_int)
    gamma = zeros(d_int)
    mu = Statistics.mean(x)
    std = Statistics.std(x)

    # Auto-covariance and variances for each blocking transformation
    for i in 1:d_int
        # n changes in length
        n = length(x)

        # Autocovariance of x
        gamma[i] = (1 / n) * sum((x[1:n-1] .- mu) .* (x[2:n] .- mu))

        # Variance of x
        s[i] = Statistics.var(x)

        # We might get a situaion where the array is not easily split in two equal sizes
        x_1 = x[1:2:end] # Extracting all numbers at odd positions
        x_2 = x[2:2:end] # Numbers at even positions
        # If length is not equal, remove highest number
        if (length(x_1) > length(x_2))
            x_1 = x_1[1:end-1]
        elseif (length(x_2) > length(x_1))
            x_2 = x_2[1:end-1]
        end
        # Blocking transformation
        x = 0.5 .* (x_1 .+ x_2)
    end

    # Test observator from theorem (chi^2-distributed)
    factor_1 = (gamma ./ s).^2
    factor_2 = 2 .^[i for i in 1:d_int]
    # Do the same length check again
    if (length(factor_1) > length(factor_2))
        factor_1 = factor_1[1:end-1]
    elseif (length(factor_2) > length(factor_1))
        factor_2 = factor_2[1:end-1]
    end

    M = (cumsum((factor_1 .* factor_2[end:-1:1])[end:-1:1]))[end:-1:1]

    # Test percentiles
    q = [6.634897,  9.210340, 11.344867, 13.276704, 15.086272,
        16.811894, 18.475307, 20.090235, 21.665994, 23.209251,
        24.724970, 26.216967, 27.688250, 29.141238, 30.577914,
        31.999927, 33.408664, 34.805306, 36.190869, 37.566235,
        38.932173, 40.289360, 41.638398, 42.979820, 44.314105,
        45.641683, 46.962942, 48.278236, 49.587884, 50.892181]

    # The actual Chi squared test - should we have stopped blocking?
    k_end = 1
    for k in 1:d_int
        k_end = k
        if (M[k] < q[k])
            break
        end
        if (k >= d)
            print("More data is needed!")
        end
    end
    
    result = s[k_end] / 2^(d-(k_end-1))

    return mu, result^0.5 * n0^0.5, std
end