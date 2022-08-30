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
    (; dims, num) = wf
    # Choosing particle to move and dimension to move it in
    new_idx = rand(rng, 1:num*dims)
    
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
    (; rng) = walker
    (; move, greens, newQF, oldQF) = walker.metro_muts
    (; dims, num) = wf
    
    # Choosing particle to move
    idx_start = rand(rng, 0:num-1) * dims + 1
    new_idx = idx_start:(idx_start + dims - 1)
    
    # Computing the quantum force for the given particle
    oldQF = QF!(oldQF, walker.positions, new_idx, wf)
    
    # The direction and length to move the chosen particle
    move .= 0.5 * metro.time_step .* oldQF
    for i in 1:dims
        move[i] += Random.randn(rng) * sqrt(metro.time_step)
        greens[i] = -move[i]
    end

    # Compute the values needed to consider accepting the move
    ratio, newQF = consider!(newQF, walker, wf, new_idx, move)

    # Evaluating the greens function
    greens .= greens .+ metro.time_step .* (oldQF .- newQF)
    greens .= greens .* (oldQF .+ newQF)
    greensFuncRatio = exp(0.5 * sum(greens))
    
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