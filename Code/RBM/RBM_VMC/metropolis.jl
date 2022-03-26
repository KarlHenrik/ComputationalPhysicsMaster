struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end

function metro_step!(particles, wf, old_amp, rng, metro::Metropolis)
    dims, num = particles.dims, particles.num
    
    # Choosing particle to move
    idx = rand(rng, 1:num * dims)
    old_pos = particles.positions[idx]
    
    # Moving chosen particle
    adj_sign = Random.rand(rng, (-1.0, 1.0))
    particles.positions[idx] += metro.step_length * adj_sign
    
    # Evaluating amplitudes
    new_amp = evaluate(particles, wf)
    wfratio = (new_amp / old_amp)^2
    
    # Accepting/Denyting new state
    if (Random.rand(rng) < wfratio)
        return true, new_amp
    else
        particles.positions[idx] = old_pos
        return false, old_amp
    end
end

struct Importance
    equil_steps::Int64
    sample_steps::Int64
    time_step::Float64
end

function metro_step!(particles, wf, old_amp, rng, metro::Importance)
    temp_vec = particles.temp_vec
    oldQF = particles.temp_vec2
    newQF = particles.temp_vec3
    dims, num = particles.dims, particles.num
    
    # Choosing particle to move
    idx_start = rand(rng, 0:num-1) * dims + 1
    idx = idx_start:idx_start + dims - 1
    old_pos = particles.positions[idx]
    QF!(oldQF, particles, idx, wf)
    
    # Moving chosen particle
    temp_vec .= 0.5 * metro.time_step .* oldQF
    temp_vec .= temp_vec .+ Random.randn(rng, dims) .* sqrt(metro.time_step)
    particles.positions[idx] .+= temp_vec
    
    # Evaluating amplitudes
    new_amp = evaluate(particles, wf)
    wfratio = (new_amp / old_amp)^2
    QF!(newQF, particles, idx, wf)
    
    # Evaluating the greens function
    temp_vec .= old_pos .- particles.positions[idx]
    temp_vec .= temp_vec .+ metro.time_step .* (oldQF .- newQF)
    temp_vec .= temp_vec .* (oldQF .+ newQF)
    greensFuncRatio = exp(0.5 * sum(temp_vec))
    
    # Accepting/Denyting new state
    if (Random.rand(rng) < greensFuncRatio * wfratio)
        return true, new_amp
    else
        particles.positions[idx] .= old_pos
        return false, old_amp
    end
end