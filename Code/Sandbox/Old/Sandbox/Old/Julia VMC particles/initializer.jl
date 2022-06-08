import Random
import LinearAlgebra as la
import StaticArrays

function initialize(num, dims, rng::Random.AbstractRNG, radius = 0)
    """
    Returns a (num, dims) Float64 Array with positions from the random uniform distribution.
    
    The `RandomDevice` `rng` is used for random number generation. The particles are at least a distance `radius` apart.
    """
    positions = zeros(dims, num)
    
    for new_particle in 1:num
        notFound = true
        attempts = 0
        local candidate
        while notFound
            candidate = (rand(rng, Float64, dims) .- 0.5) .* 2
            notFound = false
            for placed_particle in 1:new_particle - 1
                if la.norm(positions[:, placed_particle] .- candidate) < radius
                    notFound = true
                end
            end
            if attempts > 100
                error("No more room for particles. Reduce the number of particles or the hard-shell radius")
            end
            attempts += 1
        end
        positions[:, new_particle] = candidate
    end
    return positions
end