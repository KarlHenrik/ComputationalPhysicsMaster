struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end

function metro_step!(state, system, metro::Metropolis)
    # Choosing particle to move and dimension to move it in
    new_particle = rand(state.rng, 1:system.num)
    new_dim = rand(state.rng, 1:system.dims)
    
    # The direction and length to move the chosen particle
    distance = Random.rand(state.rng, (-1.0, 1.0)) * metro.step_length
    
    # Compute the values needed to consider accepting the move
    ratio_ = consider!(state, system, new_particle, new_dim, distance)

    # Accepting/Denyting new state
    if (Random.rand(state.rng) < ratio_)
        # Accept the move and update all the needed values
        accept!(state, system, new_idx, move)
        return state
    else
        return state
    end
end

struct Importance
    equil_steps::Int64
    sample_steps::Int64
    time_step::Float64
end

function metro_step!(state, system, metro::Importance)
    move = state.metro_m.move
    greens = state.metro_m.greens
    oldQF = state.metro_m.oldQF
    newQF = state.metro_m.oldQF
    
    # Choosing particle to move
    idx_start = rand(state.rng, 0:system.num-1) * system.dims + 1
    new_idx = idx_start:idx_start + system.dims - 1
    
    # Computing the quantum force for the given particle
    oldQF .= QF!(oldQF, state.positions, new_idx, system.wf)
    
    # The direction and length to move the chosen particle
    move .= 0.5 * metro.time_step .* state.metro_m.oldQF
    move .= move .+ Random.randn(state.rng, system.dims) .* sqrt(metro.time_step)
    
    # Compute the values needed to consider accepting the move
    ratio = consider!(state, system, new_idx, move)
    
    # Evaluating the greens function
    greens .= -move
    greens .= greens .+ metro.time_step .* (state.metro_m.oldQF .- state.metro_m.newQF)
    greens .= greens .* (state.metro_m.oldQF .+ state.metro_m.newQF)
    greensFuncRatio = exp(0.5 * sum(greens))
    
    # Accepting/Denyting new state
    if (Random.rand(state.rng) < greensFuncRatio * ratio)
        # Accept the move and update all the needed values
        accept!(state, system, new_idx, move)
        return state
    else
        return state
    end
end