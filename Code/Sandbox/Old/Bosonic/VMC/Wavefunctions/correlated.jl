struct Correlated{V} <: WaveFunction
    alpha::Float64
    a::Float64
    HOshape::V
    HOshape2::V
    function Correlated(alpha, a, HOshape)
        HOshape = sa.SVector{size(HOshape,1), Float64}(HOshape)
        return new{typeof(HOshape)}(alpha, a, HOshape, HOshape.^2)
    end
end


function ratio(particles, p1, old_pos, wf::Correlated)::Float64
    temp_vec = particles.temp_vec
    new_pos = particles.positions[p1]
    
    temp_vec .= old_pos.^2
    temp_vec .= temp_vec .- new_pos.^2
    temp_vec .= temp_vec .* wf.HOshape
    new_term = exp(wf.alpha * sum(temp_vec))
    old_term = 1
    
    for p2 in 1:particles.num
        if p2 != p1
            p2_pos = particles.positions[p2]
            temp_vec .= p2_pos .- new_pos
            dist = la.norm(temp_vec)
            if dist < wf.a
                return 0.0
            end
            new_term *= 1.0 - wf.a / dist
            
            temp_vec .= p2_pos .- old_pos
            old_term *= 1.0 - wf.a / la.norm(temp_vec)
        end
    end
    
    return new_term / old_term
end

function kinetic(particles, wf::Correlated)::Float64
    dims, num = particles.dims, particles.num
    temp_vec = zero(particles.temp_vec)
    dblDer = 0
    
    # uncorrelated part
    for pos in particles.positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* wf.HOshape2
    dblDer += 4 * wf.alpha^2 * sum(temp_vec)
    dblDer += -2 * num * wf.alpha * sum(wf.HOshape)
    
    # vector sum and final sum
    for p1 in 1:num
        p1_pos = particles.positions[p1]
        vecSum = zero(temp_vec) # How does this not allocate!?
        for p2 in 1:num
            if p1 != p2
                temp_vec .= particles.positions[p2] .- p1_pos
                r = la.norm(temp_vec)
                fac = wf.a / (r^3 - wf.a * r^2)
                
                vecSum .+= temp_vec .* fac
                dblDer += (wf.a^2 - 2 * wf.a * r) / (r^2 - wf.a * r)^2 + fac # the final sum
            end
        end
        # adding the vector products
        temp_vec .= wf.HOshape .* p1_pos
        temp_vec .= temp_vec .* vecSum
        vecSum .= vecSum.^2 # temporarily storing the squared elements
        dblDer += sum(vecSum) + 4 * wf.alpha * sum(temp_vec)
    end
    
    return -0.5 * dblDer #-0.5 to go from double derivative to kinetic energy
end

function QF(particles, p1, wf::Correlated)
    temp_vec = zero(particles.temp_vec)
    num = particles.num
    new_pos = particles.positions[p1]
    
    temp_vec .= new_pos .* wf.HOshape
    vecSum = -4 * wf.alpha * temp_vec
    for p2 in 1:num
        if p2 != p1
            temp_vec .= particles.positions[p2] .- new_pos
            r = la.norm(temp_vec)
            fac = 2 * wf.a / (r^3 - wf.a * r^2)
            vecSum .+= temp_vec .* fac
        end
    end
    
    return vecSum
end

function paramDer(particles, wf::Correlated)::Float64
    temp_vec = zero(particles.temp_vec) #maybe there is a faster way? why does this not allocate??
    for pos in particles.positions
        temp_vec .+= pos.^2
    end
    temp_vec .= temp_vec .* wf.HOshape
    return -sum(temp_vec) / particles.num
end

function applyGradient(wf::Correlated, grad)
    return Correlated(wf.alpha - grad, wf.a, wf.HOshape)
end;