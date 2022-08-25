struct Correlated <: WaveFunction
    α::Float64
    a::Float64
    dims::Int64
    num::Int64
    HOshape::Vector{Float64}
    HOshape2::Vector{Float64}
    temp_vec::Vector{Float64}
    function Correlated(dims, num; α, a, HOshape=ones(dims))
        @assert dims == length(HOshape)
        return new(α, a, dims, num, HOshape, HOshape.^2, HOshape.*0)
    end
end
private_wf(wf::Correlated) = Correlated(wf.dims, wf.num, α=wf.α, a=wf.a, HOshape=wf.HOshape)

function ratio_direct(wf::Correlated, positions, new_idx::Int64, old_pos)
    (; alpha, α, a, HOshape, temp_vec, num) = wf
    new_pos = positions[new_idx]
    
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- new_pos.^2
    temp_vec .= temp_vec .* HOshape
    new_term = exp(α * sum(temp_vec))
    old_term = 1
    
    for p2 in 1:num
        if p2 != p1
            p2_pos = positions[p2]
            temp_vec .= p2_pos .- new_pos
            dist = la.norm(temp_vec)
            if dist < a
                return 0.0
            end
            new_term *= 1.0 - a / dist
            
            temp_vec .= p2_pos .- old_pos
            old_term *= 1.0 - a / la.norm(temp_vec)
        end
    end
    
    return new_term / old_term
end

function kinetic(particles, wf::Correlated)::Float64
    (; dims, num, temp_vec, HOshape, α, a) = wf
    dblDer = 0
    
    # uncorrelated part
    for pos in positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* HOshape2
    dblDer += 4 * α^2 * sum(temp_vec)
    dblDer += -2 * num * α * sum(HOshape)
    
    # vector sum and final sum
    for p1 in 1:num
        p1_pos = positions[p1]
        vecSum = zero(temp_vec) # How does this not allocate!?
        for p2 in 1:num
            if p1 != p2
                temp_vec .= positions[p2] .- p1_pos
                r = la.norm(temp_vec)
                fac = wf.a / (r^3 - a * r^2)
                
                vecSum .+= temp_vec .* fac
                dblDer += (a^2 - 2 * a * r) / (r^2 - a * r)^2 + fac # the final sum
            end
        end
        # adding the vector products
        temp_vec .= HOshape .* p1_pos
        temp_vec .= temp_vec .* vecSum
        vecSum .= vecSum.^2 # temporarily storing the squared elements
        dblDer += sum(vecSum) + 4 * α * sum(temp_vec)
    end
    
    return -0.5 * dblDer #-0.5 to go from double derivative to kinetic energy
end

function QF!(qf, positions, idx, wf::Correlated)
    (; num, dims, a) = wf
    temp_vec = zero(wf.temp_vec)
    new_pos = positions[idx]
    p1 = (new_idx-1)÷dims+1
    
    temp_vec .= new_pos .* wf.HOshape
    qf .= -4 * wf.alpha * temp_vec
    for p2 in 1:num
        if p2 != p1
            temp_vec .= positions[(p2-1)*dims+1:(p2-1)*dims+3] .- new_pos
            r = la.norm(temp_vec)
            fac = 2 * a / (r^3 - a * r^2)
            qf .+= temp_vec .* fac
        end
    end
    
    return qf
end

function paramDer(positions, wf::Correlated)::Float64
    temp_vec = zero(wf.temp_vec) #maybe there is a faster way? why does this not allocate??
    for pos in positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* wf.HOshape
    return -sum(temp_vec) / wf.num
end

function applyGradient(wf::Correlated, grad)
    return Correlated(wf.dims, wf.num, α=wf.α - grad, a=wf.a, HOshape=wf.HOshape)
end;