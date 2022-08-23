abstract type Hamiltonian end

struct HarmonicOscillator <: Hamiltonian
    ω2::Float64
    dims::Int64
    num::Int64
    HOshape::Vector{Float64}
    function HarmonicOscillator(HOshape, num, ω)
        return new(ω^2, length(HOshape), num, HOshape)
    end
end
HarmonicOscillator(HOshape, num; ω) = HarmonicOscillator(HOshape, num, ω)
HarmonicOscillator(dims::Int64, num; ω) = HarmonicOscillator(ones(dims), num, ω)


function potential(positions, ham::HarmonicOscillator)::Float64
    """
    V = Σ 0.5 * ω^2 * r_i^2
    
    with support for elliptical wells
    """
    pot = 0.0
    for d in 1:ham.dims
        for i in d:ham.dims:length(positions)
            pot += positions[i]^2 * ham.HOshape[d]
        end
    end
    return 0.5 * ham.ω2 * pot
end

struct HORepulsion <: Hamiltonian
    ω2::Float64
    dims::Int64
    num::Int64
    HORepulsion(dims, num, ω) = new(dims, num, ω^2)
end
HORepulsion(dims, num; ω) = HORepulsion(dims, num, ω)

function potential(positions, ham::HORepulsion)::Float64
    """
    V = Σ 0.5 * ω^2 * r_i^2 + ΣΣ 1 / r_{i,j}
    """
    (;dims, num, ω2) = ham
    pot = 0.0
    for p1 in 1:num
        idx1 = (p1-1) * dims
        for d in 1:dims
            pot += 0.5 * ω2 * positions[idx1 + d]^2
        end
        
        for p2 in p1+1:num
            idx2 = (p2-1) * dims
            r12 = 0.0
            for d in 1:dims
                r12 += (positions[idx2 + d] - positions[idx1 + d])^2
            end
            pot += 1 / sqrt(r12)
        end
    end
    return pot
end