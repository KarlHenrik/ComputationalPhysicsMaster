mutable struct Walker{W<:Wf_m, M<:Metro_m, S<:Sample_m}
    const positions::Vector{Float64}
    accepted::Bool
    const rng::Random.MersenneTwister
    
    const wf_m::W
    const metro_m::M
    const sampler_m::S
end

function EquilWalker(system, metro)
    wf = system.wf
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)

    wf_m  = get_wf_m(positions, wf)
    metro_m = get_metro_m(positions, wf, metro)
    sample_m = get_no_sample()

    return Walker{typeof(wf_m), typeof(metro_m), typeof(sample_m)}(positions, rng, wf_m, metro_m, sample_m)
end

function SampledWalker(walker, system, metro, scheme)
    wf = system.wf
    rng = Random.MersenneTwister()
    positions = walker.positions

    wf_m  = get_wf_m(positions, wf)
    metro_m = get_metro_m(positions, wf, metro)
    sample_m = get_sample_m(positions, wf, scheme)

    return Walker{typeof(wf_m), typeof(metro_m), typeof(sample_m)}(positions, rng, wf_m, metro_m, sample_m)
end

end

function getPositions(rng, wf)
    # Uniform random positions in the range [-1, 1] in each dimention
    positions = (Random.rand(rng, Float64, wf.num * wf.dims) .- 0.5) .* 2
    return positions
end

function getPositions(rng, wf::Correlated)
    # Uniform random positions in the range [-1, 1] in each dimention. No particles are closer than the hard shell radius a
    dims = wf.dims
    num = wf.num
    positions = zeros(dims, num)

    for new_particle in 1:num
        tooClose = true
        attempts = 0
        local candidate
        while tooClose
            attempts += 1
            if attempts > 1000
                error("No more room for particles. Reduce the number of particles or the hard-shell radius")
            end

            candidate = (Random.rand(rng, Float64, dims) .- 0.5) .* 2
            tooClose = false

            # Checking if any of the particles preceding the new particle is too close
            for placed_particle in 1:new_particle - 1
                idx = (placed_particle - 1) * dims + 1
                if la.norm(positions[idx:idx+dims-1] .- candidate) < wf.a
                    tooClose = true
                end
            end
            # If not, exit the while loop and set the new position
        end
        idx = (new_particle - 1) * dims + 1
        positions[idx:idx+dims-1] = candidate
    end
    return vec(positions)
end

include("walker_updating.jl")