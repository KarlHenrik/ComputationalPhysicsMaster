abstract type Hamiltonian end

struct HarmonicOscillator <: Hamiltonian
    dims::Int64
    num::Int64
    ω2::Float64
    HOshape::Vector{Float64}
    function HarmonicOscillator(dims, num; ω, HOshape=ones(dims))
        @assert dims == length(HOshape)
        return new(dims, num, ω^2, HOshape)
    end
end

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
    dims::Int64
    num::Int64
    shielding2::Float64
    ω2::Float64
    HORepulsion(dims, num, shielding, ω) = new(dims, num, shielding^2, ω^2)
end
HORepulsion(dims, num; shielding, ω) = HORepulsion(dims, num, shielding, ω)

function potential(positions, ham::HORepulsion)::Float64
    """
    V = Σ 0.5 * ω^2 * r_i^2 + ΣΣ 1 / r_{i,j}
    """
    (; dims, num, shielding, ω2) = ham
    pot = 0.0
    for p1 in 1:num
        idx1 = (p1-1) * dims
        for d in 1:dims
            pot += 0.5 * ω2 * positions[idx1 + d]^2
        end
        
        for p2 in p1+1:num
            idx2 = (p2-1) * dims
            r12 = shielding2
            for d in 1:dims
                r12 += (positions[idx2 + d] - positions[idx1 + d])^2
            end
            pot += 1 / sqrt(r12)
        end
    end
    return pot
end