import LinearAlgebra as la
import StaticArrays as sa

struct SimpleGaussian{V} <: WaveFunction
    alpha::Float64
    HOshape::V
    HOshape2::V
    function SimpleGaussian(alpha, HOshape)
        HOshape = sa.SVector{size(HOshape,1), Float64}(HOshape)
        return new{typeof(HOshape)}(alpha, HOshape, HOshape.^2)
    end
end


function ratio(particles, p1, old_pos, wf::SimpleGaussian)::Float64
    """
    Old wavefunc value term: exp(-alpha * old_r2)
    New wavefunc value term: exp(-alpha * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    temp_vec = particles.temp_vec
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- particles.positions[p1].^2
    temp_vec .= temp_vec .* wf.HOshape
    return exp(wf.alpha * sum(temp_vec))
end

function kinetic(particles, wf::SimpleGaussian)::Float64
    temp_vec = zero(particles.temp_vec) #maybe there is a faster way? why does this not allocate??
    for pos in particles.positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* wf.HOshape2
    return wf.alpha * (sum(wf.HOshape) * particles.num - 2.0 * wf.alpha * sum(temp_vec))
end

function QF(particles, p1, wf::SimpleGaussian) #maybe faster to put return type here?
    particles.temp_vec .= -4 * wf.alpha .* particles.positions[p1] .* wf.HOshape
    return particles.temp_vec
end

function paramDer(particles, wf::SimpleGaussian)::Float64
    temp_vec = zero(particles.temp_vec) #maybe there is a faster way? why does this not allocate??
    for pos in particles.positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* wf.HOshape
    return -sum(temp_vec) / particles.num
end

function applyGradient(wf::SimpleGaussian, grad)
    return SimpleGaussian(wf.alpha - grad, wf.HOshape)
end