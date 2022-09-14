mutable struct Walker{S}
    const positions::Vector{Float64}
    accepted::Bool
    const rng::Random.MersenneTwister
    new_amp::Float64
    old_amp::Float64
    old_pos::Float64
    
    const samp_muts::S
end

function EquilWalker(wf, metro) # A walker for equilibrium steps
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)
    
    samp_muts = No_Muts()
    
    return Walker{typeof(samp_muts)}(positions, false, rng, 0, 0, samp_muts)
end

function SampledWalker(wf, metro, sampler, walker) # A walker using the position of a previous walker
    rng = Random.MersenneTwister()
    positions = walker.positions
    
    samp_muts = Sample_Muts(wf, sampler)
    
    return Walker{typeof(samp_muts)}(positions, false, rng, 0, 0, samp_muts)
end

function SampledWalker(wf, metro, sampler)
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)
    
    samp_muts = Sample_Muts(wf, sampler)
    
    return Walker{typeof(samp_muts)}(positions, false, rng, 0, 0, samp_muts)
end

function getPositions(rng, wf)
    # Uniform random positions in the range [-1, 1] in each dimention
    positions = (Random.rand(rng, Float64, wf.n) .- 0.5) .* 2
    return positions
end

function getPositions(rng, wf::Correlated)
    # Uniform random positions in the range [-1, 1] in each dimention. No particles are closer than the hard shell radius a
    (; n, a) = wf
    positions = zeros(n)

    for new_particle in 1:n
        tooClose = true
        attempts = 0
        local candidate
        while tooClose
            attempts += 1
            if attempts > 1000
                error("No more room for particles. Reduce the number of particles or the hard-shell radius")
            end

            candidate = (Random.rand(rng) - 0.5) * 2
            tooClose = false

            # Checking if any of the particles preceding the new particle is too close
            for placed_particle in 1:new_particle - 1
                if abs(positions[placed_particle] - candidate) < a
                    tooClose = true
                end
            end
            # If not, exit the while loop and set the new position
        end
        positions[new_particle] = candidate
    end
    return positions
end