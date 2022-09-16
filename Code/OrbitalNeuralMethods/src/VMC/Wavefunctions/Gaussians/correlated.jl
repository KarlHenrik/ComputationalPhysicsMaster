struct Correlated <: Gaussian
    α::Float64
    a::Float64
    n::Int64
    
    function Correlated(n; α, a)
        return new(α, a, n)
    end
end
private_wf(wf::Correlated, positions) = Correlated(wf.n, α=wf.α, a=wf.a)

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
    (; n, α, a) = wf
    
    dblDer = 0.0
    
    kin_sum = 0.0
    for pos in positions
        kin_sum += pos^2
    end
    dblDer += -2 * α * (n - 2 * α * kin_sum)
    
    # vector sum and final sum
    @inbounds for p1 in 1:n
        p1_pos = positions[p1]
        vecSum = 0.0
        for p2 in 1:n
            if p1 != p2
                r = abs(positions[p2] - p1_pos)
                fac = a / (r^3 - a * r^2)
                
                vecSum += temp_vec .* fac
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

function applyGradient(wf::Correlated, grad)
    return Correlated(wf.n, α=wf.α - grad, a=wf.a)
end;