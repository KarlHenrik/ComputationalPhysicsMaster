using Random
include("Wavefunctions/wavefunctions.jl")

struct Metropolis
    equil_steps::Int64
    sample_steps::Int64
    step_length::Float64
end

function metro_step!(particles, wf, rng, metro::Metropolis)::Bool
    dims, num = particles.dims, particles.num
    p1 = rand(rng, 1:num)
    old_pos = copy(particles.positions[p1])

    adj_dir = rand(rng, 1:dims)
    adj_sign = rand(rng, (-1, 1))
    particles.positions[p1][adj_dir] += metro.step_length * adj_sign
    
    wfratio = ratio(particles, p1, old_pos, wf)
    if (rand(rng) < wfratio^2)
        return true
    else
        particles.positions[p1] .= old_pos
        return false
    end
end

struct Importance
    equil_steps::Int64
    sample_steps::Int64
    time_step::Float64
end

function metro_step!(particles, wf, rng, metro::Importance)::Bool
    temp_vec = particles.temp_vec
    dims, num = particles.dims, particles.num
    p1 = rand(rng, 1:num)
    old_pos = copy(particles.positions[p1])
    oldQF = QF(particles, p1, wf)

    temp_vec .= 0.5 * metro.time_step .* oldQF
    temp_vec .= temp_vec .+ randn(rng, dims) .* sqrt(metro.time_step)
    particles.positions[p1] .+= temp_vec
    
    wfratio = ratio(particles, p1, old_pos, wf)
    newQF = QF(particles, p1, wf)
    
    temp_vec .= old_pos .- particles.positions[p1]
    temp_vec .= temp_vec .+ metro.time_step .* (oldQF .- newQF)
    temp_vec .= temp_vec .* (oldQF .+ newQF)
    greensFuncRatio = exp(0.5 * sum(temp_vec))
    
    if (rand(rng) < greensFuncRatio * wfratio^2)
        return true
    else
        particles.positions[p1] .= old_pos
        return false
    end
end

