function potential(positions, ham::HarmonicOscillator)::Float64
    """
    V = Σ 0.5 * ω^2 * r_i^2
    
    with support for elliptical wells
    """
    pot = 0.0
    for pos in positions
        pot += pos^2
    end
    
    return 0.5 * ham.ω2 * pot
end

function potential(positions, ham::HOCoulomb)::Float64
    """
    V = Σ 0.5 * ω^2 * r_i^2 + ΣΣ 1 / r_{i,j}
    """
    (; shielding2, ω2) = ham
    pot = 0.0
    for p1 in eachindex(positions)
        pot += 0.5 * ω2 * positions[p1]^2
        
        for p2 in p1+1:length(positions)
            r12 = shielding2 + (positions[p2] - positions[p1])^2

            pot += 1 / sqrt(r12)
        end
    end
    return pot
end