struct Correlated <: WaveFunction
    α::Float64
    a::Float64
    dims::Int64
    num::Int64
    HOshape::Vector{Float64}
    HOshape2::Vector{Float64}
    p1_pos::Vector{Float64}
    p2_pos::Vector{Float64}
    temp_vec::Vector{Float64}
    function Correlated(dims, num; α, a, HOshape=ones(dims))
        @assert dims == length(HOshape)
        return new(α, a, dims, num, HOshape, HOshape.^2, HOshape.*0, HOshape.*0, HOshape.*0)
    end
end
private_wf(wf::Correlated, positions) = Correlated(wf.dims, wf.num, α=wf.α, a=wf.a, HOshape=copy(wf.HOshape))

# Metropolis version
function ratio_direct(wf::Correlated, positions, new_idx::Int64, old_pos_i)
    (; α, a, HOshape, temp_vec, num, dims, p1_pos, p2_pos) = wf
    p1 = new_idx - ((new_idx+(dims-1))%dims)
    old_dim = (new_idx+dims-1)%dims+1

    for i in 1:dims
        p1_pos[i] = positions[(p1:p1+dims-1)[i]]
        temp_vec[i] = 0.0
    end
    temp_vec[old_dim] = old_pos_i^2 - p1_pos[old_dim]^2
    temp_vec .= temp_vec .* HOshape
    
    new_term = exp(α * sum(temp_vec))
    old_term = 1.0
    
    for p2 in 1:dims:num*dims
        if p2 != p1
            for (j, p2_i) in enumerate(p2:p2+dims-1)
                p2_pos[j] = positions[p2_i]
                temp_vec[j] = p2_pos[j] - p1_pos[j]
            end
            
            dist = la.norm(temp_vec)
            if dist < a
                return 0.0
            end
            new_term *= 1.0 - a / dist
            
            temp_vec .= p2_pos .- p1_pos
            temp_vec[old_dim] = p2_pos[old_dim] - old_pos_i
            old_term *= 1.0 - a / la.norm(temp_vec)
        end
    end
    
    return new_term / old_term
end

# Importance version
function ratio_direct(wf::Correlated, positions, new_idx::UnitRange{Int64}, old_pos)
    (; α, a, HOshape, temp_vec, num, dims, p1_pos, p2_pos) = wf
    p1 = new_idx[1]
    for i in 1:dims
        p1_pos[i] = positions[new_idx[i]]
    end
    
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- p1_pos.^2
    temp_vec .= temp_vec .* HOshape
    
    new_term = exp(α * sum(temp_vec))
    old_term = 1.0
    
    for p2 in 1:dims:num*dims
        if p2 != p1
            for (j, p2_i) in enumerate(p2:p2+dims-1)
                p2_pos[j] = positions[p2_i]
                temp_vec[j] = p2_pos[j] - p1_pos[j]
            end
            
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
    (; dims, num, temp_vec, HOshape, HOshape2, α, a, p1_pos) = wf
    vecSum = wf.p2_pos
    
    dblDer = 0.0
    
    kin_sum = 0.0
    @inbounds for d in 1:dims
        for i in d:dims:length(positions)
            kin_sum += positions[i]^2 * HOshape2[d]
        end
    end
    dblDer += -2 * α * (sum(wf.HOshape) * wf.num - 2 * α * kin_sum)
    
    # vector sum and final sum
    @inbounds for p1 in 1:dims:num*dims
        for i in 1:dims
            p1_pos[i] = positions[p1 + i - 1]
            vecSum[i] = 0.0
        end
        for p2 in 1:dims:num*dims
            if p1 != p2
                for i in 1:dims
                    temp_vec[i] = positions[p2 + i - 1] .- p1_pos[i]
                end
                
                r = la.norm(temp_vec)
                fac = a / (r^3 - a * r^2)
                
                vecSum .+= temp_vec .* fac
                dblDer += (a^2 - 2 * a * r) / (r^2 - a * r)^2 + 2 * fac # the final sum
            end
        end
        # adding the vector products
        temp_vec .= HOshape .* p1_pos .* vecSum
        vecSum .= vecSum.^2 # temporarily storing the squared elements
        dblDer += sum(vecSum) + 4 * α * sum(temp_vec)
    end
    
    return -0.5 * dblDer #-0.5 to go from double derivative to kinetic energy
end

function QF!(qf, positions, idx, wf::Correlated)
    (; num, dims, a, α, temp_vec, p1_pos, HOshape) = wf
    for i in 1:dims
        p1_pos[i] = positions[idx[i]]
        qf[i] = -4 * α * p1_pos[i] * HOshape[i]
    end
    p1 = idx[1]
    
    for p2 in 1:dims:num*dims
        if p2 != p1
            for i in 1:dims
                temp_vec[i] = positions[p2+i-1] - p1_pos[i]
            end
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