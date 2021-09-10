abstract type Hamiltonian end

struct HarmonicOscillator<:Hamiltonian
    omega2::Float64
    HOshape::Array{Float64}
    HarmonicOscillator(omega, HOshape) = new(omega^2, HOshape)
end

function potential(particles, ham::HarmonicOscillator)::Float64
    particles.pos_temp .= particles.positions.^2
    particles.pos_temp .= particles.pos_temp .* ham.HOshape
    return 0.5 * ham.omega2 * sum(particles.pos_temp)
end