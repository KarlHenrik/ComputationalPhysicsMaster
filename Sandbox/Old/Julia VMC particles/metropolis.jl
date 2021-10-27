using Random
include("wavefunctions.jl")

function metro_step!(particles, wf, rng, step_length)::Bool
    p1 = rand(rng, 1:particles.num)
    old_pos = particles.positions[:, p1]

    adj_dir = rand(rng, 1:particles.dims)
    adj_sign = rand(rng, (-1, 1))
    new_pos = copy(old_pos)
    new_pos[adj_dir] += step_length * adj_sign

    particles.positions[:, p1] .= new_pos
    particles.pos2[:, p1] .= new_pos.^2
    wfratio = ratio(particles, p1, old_pos, wf)
    if (rand(rng) < wfratio^2)
        return true
    else
        particles.positions[:, p1] .= old_pos
        particles.pos2[:, p1] .= old_pos.^2
        return false
    end
    return
end

#=
function importance_step(positions, wf, rng, num, dims, step_length):
    p1 = tf.math.argmax(chosen_particle_mask)[0]
    oldQF = wf.QF(old_positions, p1)

    move = 0.5 * oldQF * step_length + tf.random.normal(shape = (dims,), dtype = tf.float64) * tf.math.sqrt(step_length)
    new_positions = old_positions + chosen_particle_mask * move

    ratio = wf.ratio(new_positions, p1, old_positions)
    newQF = wf.QF(new_positions, p1)

    greensFuncRatio = tf.math.reduce_sum( (oldQF + newQF) * (old_positions[p1] - new_positions[p1] + step_length * (oldQF - newQF)) )
    greensFuncRatio = tf.math.exp(0.5 * greensFuncRatio)
    if (tf.random.uniform(shape = (), dtype = tf.float64) < greensFuncRatio * tf.math.pow(ratio, 2)):
        return new_positions
    else:
        return old_positions
=#


