struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end
function Metropolis(;equils, samples, step)
    return Metropolis(equils, samples, step)
end

function metro_step!(walker, wf, metro::Metropolis)
    (; rng, positions) = walker
    (; n) = wf
    
    new_idx = rand(rng, 1:n) # Choosing particle to move and dimension to move it in
    move = Random.rand(rng, (-1.0, 1.0)) * metro.step_length # The direction and length to move the chosen particle
    
    old_pos = positions[new_idx]
    positions[new_idx] += move

    # Compute the values needed to consider accepting the move
    ratio = consider!(wf, positions, new_idx, old_pos)
    
    # Accepting/Denyting new walker
    if (Random.rand(rng) < ratio^2)
        # Update values to be sampled (And saved NN QF values)
        accept!(wf, new_idx, positions[new_idx])
        return true
    else
        # Revert positions (And selected Slater matrix columns)
        positions[new_idx] = old_pos
        return false
    end
end

struct Importance
    equil_steps::Int64
    sample_steps::Int64
    time_step::Float64
end
function Importance(;equils, samples, step)
    return Importance(equils, samples, step)
end

function metro_step!(walker, wf, metro::Importance)
    (; rng, positions) = walker
    
    # Generating a move
    new_idx = rand(rng, 1:wf.n) # Choosing particle to move
    oldQF = QF(positions, new_idx, wf) # Computing the quantum force for the given particle
    move = 0.5 * metro.time_step * oldQF
    move += Random.randn(rng) * sqrt(metro.time_step) # The direction and length to move the chosen particle

    # Moving the particle
    old_pos = positions[new_idx]
    positions[new_idx] += move
    ratio, newQF = consider_qf!(wf, positions, new_idx, old_pos) # Compute the values needed to consider accepting the move

    # Evaluating the greens function
    greens = -move
    greens += metro.time_step * (oldQF - newQF)

    greens = greens * (oldQF + newQF)
    greensFuncRatio = exp(0.5 * greens)
    
    # Accepting/Denyting new position
    if (Random.rand(rng) < greensFuncRatio * ratio^2)
        # Update values to be sampled (And saved NN QF values)
        wf = accept!(wf, new_idx, positions[new_idx])
        return true
    else
        # Revert positions (And selected Slater matrix columns)
        positions[new_idx] = old_pos
        return false
    end
end