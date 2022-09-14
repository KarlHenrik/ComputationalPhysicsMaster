struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end
function Metropolis(;equils, samples, step)
    return Metropolis(equils, samples, step)
end

function metro_step!(walker, wf, metro::Metropolis, ham)
    (; rng) = walker
    (; n) = wf
    # Choosing particle to move and dimension to move it in
    new_idx = rand(rng, 1:n)
    
    # The direction and length to move the chosen particle
    move = Random.rand(rng, (-1.0, 1.0)) * metro.step_length
    
    # Compute the values needed to consider accepting the move
    ratio = consider!(walker, wf, new_idx, move)
    
    # Accepting/Denyting new walker
    if (Random.rand(rng) < ratio^2)
        # Update values to be sampled (And saved NN QF values)
        accept!(walker, wf, new_idx, ham)
        return walker
    else
        # Revert positions (And selected Slater matrix columns)
        deny!(walker, wf, new_idx)
        return walker
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

function metro_step!(walker, wf, metro::Importance, ham)
    (; rng, positions) = walker
    (; n) = wf
    
    # Choosing particle to move
    new_idx = rand(rng, 1:n)
    
    # Computing the quantum force for the given particle
    oldQF = QF(positions, new_idx, wf)
    
    # The direction and length to move the chosen particle
    move = 0.5 * metro.time_step * oldQF
    move += Random.randn(rng) * sqrt(metro.time_step)
    greens = -move

    # Compute the values needed to consider accepting the move
    ratio, newQF = consider!(walker, wf, new_idx, move)

    # Evaluating the greens function
    greens = greens + metro.time_step * (oldQF - newQF)
    greens = greens * (oldQF + newQF)
    greensFuncRatio = exp(0.5 * greens)
    
    # Accepting/Denyting new walker
    if (Random.rand(rng) < greensFuncRatio * ratio^2)
        # Update values to be sampled (And saved NN QF values)
        accept!(walker, wf, new_idx, ham)
        return walker
    else
        # Revert positions (And selected Slater matrix columns)
        deny!(walker, wf, new_idx)
        return walker
    end
end