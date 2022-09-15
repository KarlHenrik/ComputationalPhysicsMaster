struct Walker
    positions::Vector{Float64}
    rng::Random.MersenneTwister
end

# TODO set up distance matrix support for coulomb interacting systems
function Walker(wf)
    rng = Random.MersenneTwister()
    positions = getPositions(rng, wf)
    
    return Walker(positions, rng)
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