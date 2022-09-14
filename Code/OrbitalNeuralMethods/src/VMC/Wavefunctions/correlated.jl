struct Correlated <: WaveFunction
    α::Float64
    a::Float64
    n::Int64
    
    p1_pos::Vector{Float64}
    p2_pos::Vector{Float64}
    temp_vec::Vector{Float64}
    function Correlated(n; α, a)
        @assert dims == length(HOshape)
        return new(α, a, dims, n, HOshape, HOshape.^2, HOshape.*0, HOshape.*0, HOshape.*0)
    end
end
private_wf(wf::Correlated, positions) = Correlated(wf.dims, wf.n, α=wf.α, a=wf.a, HOshape=copy(wf.HOshape))

# Metropolis version
function ratio_direct(wf::Correlated, positions, p1::Int64, old_pos::Float64)
    (; n, α, a) = wf

    p1_pos = positions[p1]

    new_term = exp(α * (old_pos^2 - p1_pos^2))
    old_term = 1.0
    
    for p2 in 1:n
        if p2 != p1
            p2_pos = positions[p2]
            
            dist = abs(p2_pos - p1_pos)
            if dist < a
                return 0.0
            end
            new_term *= 1.0 - a / dist
            old_term *= 1.0 - a / abs(p2_pos - old_pos)
        end
    end
    
    return new_term / old_term
end

function kinetic(positions, wf::Correlated)::Float64
    (; n, temp_vec, HOshape, HOshape2, α, a, p1_pos) = wf
    
    dblDer = 0.0
    
    kin_sum = 0.0
    for pos in positions
        kin_sum += pos^2
    end
    dblDer += -2 * α * (n - 2 * α * kin_sum)
    
    # vector sum and final sum
    @inbounds for p1 in 1:n
        p1_pos = positions[p1]
        for p2 in 1:n
            if p1 != p2
                r = abs(positions[p2] - p1_pos)
                fac = a / (r^3 - a * r^2)
                
                vecSum .+= temp_vec .* fac
                dblDer += (a^2 - 2 * a * r) / (r^2 - a * r)^2 + 2 * fac # the final sum
            end
        end
        # adding the vector products
        temp_vec = p1_pos * vecSum
        vecSum = vecSum^2 # temporarily storing the squared elements
        dblDer += vecSum + 4 * α * temp_vec
    end
    
    return -0.5 * dblDer #-0.5 to go from double derivative to kinetic energy
end

function QF(positions, p1::Int64, wf::Correlated)
    (; n, a, α) = wf
    p1_pos = positions[p1]
    qf = -4 * α * p1_pos
    
    for p2 in 1:n
        if p2 != p1
            
            temp_vec = positions[p2] - p1_pos
            r = abs(temp_vec)
            fac = 2 * a / (r^3 - a * r^2)
            qf += temp_vec * fac
        end
    end
    
    return qf
end

function paramDer!(samp_muts, positions, wf::Correlated)
    pos_sum = 0.0
    for pos in positions
        pos_sum += pos^2
    end

    samp_muts.paramDer = -pos_sum / wf.n
    return samp_muts
end

function paramDerHolder(wf::Correlated)
    return 0.0
end

function applyGradient(wf::Correlated, grad)
    return Correlated(wf.n, α=wf.α - grad, a=wf.a)
end;