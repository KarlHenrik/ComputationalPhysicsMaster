abstract type WaveFunction end

struct SimpleGaussian <: WaveFunction
    alpha::Float64
    HOshape::Array{Float64}
end


function ratio(particles, p1, old_pos, wf::SimpleGaussian)::Float64
    """
    Old wavefunc value term: exp(-alpha * old_r2)
    New wavefunc value term: exp(-alpha * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    #r2 = sum(particles.positions[p1].^2 .* wf.HOshape)
    #old_r2 = sum(old_pos.^2 .* wf.HOshape)
    #return exp.(wf.alpha * (old_r2 - r2))
    old_new_ratio = sum((old_pos.^2 .- particles.positions[:, p1].^2) .* wf.HOshape)
    return exp(wf.alpha * old_new_ratio)
end

function kinetic(particles, wf::SimpleGaussian)::Float64
    particles.pos_temp .= particles.pos2 .* wf.HOshape
    particles.pos_temp .= wf.HOshape .* (1 .- 2 * wf.alpha .* particles.pos_temp)
    return wf.alpha * sum(particles.pos_temp)
end

function QF(particles, p1, wf::SimpleGaussian)::Array{Float64}
    return -4 * wf.alpha .* particles.positions[:, p1] .* wf.HOshape
end

function paramDer(particles, wf::SimpleGaussian)::Float64
    particles.pos_temp .= particles.pos2 .* wf.HOshape
    return -sum(particles.pos_temp) / particles.num
end

struct Correlated <: WaveFunction
    alpha::Float64
    a::Float64
    HOshape::Array{Float64}
end

function ratio(particles, p1, old_pos, wf::Correlated)::Float64
    old_pos.^2
end

function kinetic(particles, wf::Correlated)::Float64
    
end

function QF(particles, p1, wf::Correlated)::Array{Float64}
    
end

function paramDer(particles, wf::Correlated)::Float64
    
end