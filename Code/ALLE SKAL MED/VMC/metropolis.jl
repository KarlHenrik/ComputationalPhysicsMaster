struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end

function metro_step!(walker, system, metro::Metropolis)
    # Choosing particle to move and dimension to move it in
    new_particle = rand(walker.rng, 1:system.num)
    new_dim = rand(walker.rng, 1:system.dims)
    
    # The direction and length to move the chosen particle
    distance = Random.rand(walker.rng, (-1.0, 1.0)) * metro.step_length
    
    # Compute the values needed to consider accepting the move
    ratio_ = consider!(walker, system, new_particle, new_dim, distance)

    # Accepting/Denyting new walker
    if (Random.rand(walker.rng) < ratio_)
        # Accept the move and update all the needed values
        accept!(walker, system, new_idx, move)
        return walker
    else
        return walker
    end
end

struct Importance
    equil_steps::Int64
    sample_steps::Int64
    time_step::Float64
end

function metro_step!(walker, system, metro::Importance)
    (; move, greens, oldQF, newQF) = walker.metro_m
    
    # Choosing particle to move
    idx_start = rand(walker.rng, 0:system.num-1) * system.dims + 1
    new_idx = idx_start:idx_start + system.dims - 1
    
    # Computing the quantum force for the given particle
    oldQF .= QF!(oldQF, walker.positions, new_idx, system.wf)
    
    # The direction and length to move the chosen particle
    move .= 0.5 * metro.time_step .* walker.metro_m.oldQF
    move .= move .+ Random.randn(walker.rng, system.dims) .* sqrt(metro.time_step)
    
    # Compute the values needed to consider accepting the move
    ratio = consider!(walker, system, new_idx, move)
    
    # Evaluating the greens function
    greens .= -move
    greens .= greens .+ metro.time_step .* (oldQF .- newQF)
    greens .= greens .* (oldQF .+ newQF)
    greensFuncRatio = exp(0.5 * sum(greens))
    
    # Accepting/Denyting new walker
    if (Random.rand(walker.rng) < greensFuncRatio * ratio)
        # Accept the move and update all the needed values
        accept!(walker, system, new_idx, move)
        return walker
    else
        return walker
    end
end