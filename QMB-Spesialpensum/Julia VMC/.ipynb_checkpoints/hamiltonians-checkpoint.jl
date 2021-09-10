abstract type Hamiltonian end

struct HarmonicOscillator<:Hamiltonian
    omega2::Float64
    HOshape::Array{Float64}
    HarmonicOscillator(omega, HOshape) = new(omega^2, HOshape)
end

function potential(positions, ham::HarmonicOscillator, temp)::Float64
    temp .= positions.^2
    temp .= temp .* ham.HOshape
    return 0.5 * ham.omega2 * sum(temp)
end