abstract type Hamiltonian end

struct HarmonicOscillator <: Hamiltonian
    omega2::Float64
    HOshape::Vector{Float64}
    function HarmonicOscillator(omega, HOshape)
        return new(omega^2, HOshape)
    end
end

function potential(particles, ham::HarmonicOscillator)::Float64
    pot = 0.0
    for d in 1:particles.dims
        for i in d:particles.dims:length(particles.positions)
            pot += particles.positions[i]^2 * ham.HOshape[d]
        end
    end
    return 0.5 * ham.omega2 * pot
end

struct HORepulsion <: Hamiltonian
    omega2::Float64
    HORepulsion(omega) = new(omega^2)
end

function potential(particles, ham::HORepulsion)::Float64
    pot = 0.0
    for p1 in 1:particles.num
        idx1 = (p1-1) * particles.dims
        for d in 1:particles.dims
            pot += 0.5 * ham.omega2 * particles.positions[idx1 + d]^2
        end
        
        for p2 in p1+1:particles.num
            idx2 = (p2-1) * particles.dims
            r12 = 0.0
            for d in 1:particles.dims
                r12 += (particles.positions[idx2 + d] - particles.positions[idx1 + d])^2
            end
            pot += 1 / sqrt(r12)
        end
    end
    return pot
end