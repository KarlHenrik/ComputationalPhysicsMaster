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

function ratio_direct(wf::Correlated, positions, new_idx::Int64, old_pos_i)
    (; α, a, HOshape, temp_vec, num, dims) = wf
    p1 = new_idx - ((new_idx+(dims-1))%dims)
    
    new_pos = positions[p1:p1+dims-1]
    old_pos = copy(new_pos)
    old_pos[(new_idx-1)%dims+1] = old_pos_i
    
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- new_pos.^2
    temp_vec .= temp_vec .* HOshape
    new_term = exp(α * sum(temp_vec))
    old_term = 1.0
    
    for p2 in 1:dims:num*dims
        if p2 != p1
            p2_pos = positions[p2:p2+dims-1]
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

function ratio_direct(wf::Correlated, positions, new_idx::UnitRange{Int64}, old_pos)
    (; α, a, HOshape, temp_vec, num, dims) = wf
    p1 = new_idx[1]
    
    new_pos = positions[new_idx]
    
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- new_pos.^2
    temp_vec .= temp_vec .* HOshape
    new_term = exp(α * sum(temp_vec))
    old_term = 1.0
    
    for p2 in 1:dims:num*dims
        if p2 != p1
            p2_pos = positions[p2:p2+dims-1]
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


function kinetic(positions, wf::Correlated)::Float64
    (; dims, num, temp_vec, HOshape, HOshape2, α, a) = wf
    
    dblDer = 0.0
    
    kin_sum = 0.0
    for d in 1:dims
        for i in d:dims:length(positions)
            kin_sum += positions[i]^2 * wf.HOshape2[d]
        end
    end
    dblDer += -2 * α * (sum(wf.HOshape) * wf.num - 2 * α * kin_sum)

    # vector sum and final sum
    for p1 in 1:dims:num*dims
        p1_pos = positions[p1:p1+dims-1]
        vecSum = zero(temp_vec) # How does this not allocate!?
        for p2 in 1:dims:num*dims
            if p1 != p2
                temp_vec .= positions[p2:p2+dims-1] .- p1_pos
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
    (; num, dims, a, α, temp_vec) = wf
    temp_vec = zero(temp_vec)
    new_pos = positions[idx]
    p1 = idx[1]
    
    temp_vec .= new_pos .* wf.HOshape
    qf .= -4 * α * temp_vec
    for p2 in 1:dims:num*dims
        if p2 != p1
            temp_vec .= positions[p2:p2+dims-1] .- new_pos
            r = la.norm(temp_vec)
            fac = 2 * a / (r^3 - a * r^2)
            qf .+= temp_vec .* fac
        end
    end
    
    return qf
end

function paramDer!(samp_muts, positions, wf::Correlated)
    pos_sum = 0.0
    for d in 1:wf.dims
        for i in d:wf.dims:length(positions)
            pos_sum += positions[i]^2 * wf.HOshape[d]
        end
    end
    samp_muts.paramDer = -pos_sum / wf.num
    return samp_muts
end

function paramDerHolder(wf::Correlated)
    return 0.0
end

function applyGradient(wf::Correlated, grad)
    return Correlated(wf.dims, wf.num, α=wf.α - grad, a=wf.a, HOshape=wf.HOshape)
end;