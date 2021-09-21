abstract type Hamiltonian end

struct HarmonicOscillator{V} <: Hamiltonian
    omega2::Float64
    HOshape::V
    function HarmonicOscillator(omega, HOshape)
        HOshape = sa.SVector{size(HOshape,1), Float64}(HOshape)
        return new{typeof(HOshape)}(omega^2, HOshape)
    end
end

function potential(particles, ham::HarmonicOscillator)::Float64
    temp_vec = zero(particles.temp_vec) #maybe there is a faster way? why does this not allocate??
    for pos in particles.positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* ham.HOshape
    return 0.5 * ham.omega2 * sum(temp_vec)
end