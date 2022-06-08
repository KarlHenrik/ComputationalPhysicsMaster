struct HOBasis <: SpatialBasis
    l::Int64  # number of basis functions
    ω::Float64  # strength of harmonic oscillator potential
    hermites::Vector{Float64}
    
    function HOBasis(l, ω)
        hermites = zeros(l)
        hermites[1] = 1.0
        
        return new(l, ω, hermites)
    end
end

function evaluate(x, ho::HOBasis)
    (; ω, hermites) = ho
    hos = zero(hermites)
    
    x = √ω * x

    hermites[1] = 1.0
    hermites[2] = 2x

    ho_fac = (ω / π)^0.25 * exp(-x^2 / 2)
    hos[1] = ho_fac * hermites[1]
    ho_fac *= 1 / √2
    hos[2] = ho_fac * hermites[2]

    @inbounds for n in 3:length(hos)
        hermites[n] = 2x * hermites[n-1] - 2(n - 2) * hermites[n-2]

        ho_fac *= 1 / sqrt( 2(n - 1) )
        hos[n] = ho_fac * hermites[n]
    end
    
    return hos
end

function evaluate!(hos, x, ho::HOBasis)
    (; ω, hermites) = ho
    
    x = √ω * x

    hermites[1] = 1.0
    hermites[2] = 2x

    ho_fac = (ω / π)^0.25 * exp(-x^2 / 2)
    hos[1] = ho_fac * hermites[1]
    ho_fac *= 1 / √2
    hos[2] = ho_fac * hermites[2]

    @inbounds for n in 3:length(hos)
        hermites[n] = 2x * hermites[n-1] - 2(n - 2) * hermites[n-2]

        ho_fac *= 1 / sqrt( 2(n - 1) )
        hos[n] = ho_fac * hermites[n]
    end
    
    return hos
end

function spatial(ho::HOBasis, grid)
    n = length(grid)
    l = ho.l
    hos = zeros(l)
    res = [zeros(n) for i in 1:l]
    
    for i in 1:n
        evaluate!(hos, grid[i], ho)
        for j in 1:l
            @inbounds res[j][i] = hos[j]
        end
    end
    
    return res
end