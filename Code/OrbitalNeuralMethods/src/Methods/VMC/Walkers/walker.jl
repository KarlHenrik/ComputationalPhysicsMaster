mutable struct Walker{S, Q}
    const positions::Vector{Float64}
    accepted::Bool
    const rng::Random.MersenneTwister
    new_amp::Float64
    old_amp::Float64
    
    const samp_muts::S
    const qf_muts::Q
end

function Walker(wf, metro) # A walker for equilibrium steps
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)
    
    samp_muts = No_Muts()
    qf_muts = QF_Muts(wf, metro)
    
    return Walker{typeof(e_muts), typeof(qf_muts)}(positions, false, rng, 0, 0, e_muts, qf_muts)
end

function Walker(wf, metro, sampler, walker) # A walker using the position of a previous walker
    rng = Random.MersenneTwister()
    positions = walker.positions
    
    samp_muts = Sample_Muts(wf, sampler)
    qf_muts = QF_Muts(wf, metro)
    
    return Walker{typeof(e_muts), typeof(qf_muts)}(positions, false, rng, 0, 0, e_muts, qf_muts)
end

function Walker(wf, metro, sampler)
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)
    
    samp_muts = Sample_Muts(wf, sampler)
    qf_muts = QF_Muts(wf, metro)
    
    return Walker{typeof(e_muts), typeof(qf_muts)}(positions, false, rng, 0, 0, e_muts, qf_muts)
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