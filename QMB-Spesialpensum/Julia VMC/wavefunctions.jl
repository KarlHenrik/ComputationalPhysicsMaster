abstract type WaveFunction end

struct SimpleGaussian <: WaveFunction
    alpha::Float64
    HOshape::Array{Float64}
end


function ratio(positions, p1, old_pos, wf::SimpleGaussian)::Float64
    """
    Old wavefunc value term: exp(-alpha * old_r2)
    New wavefunc value term: exp(-alpha * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    #r2 = sum(particles.positions[p1].^2 .* wf.HOshape)
    #old_r2 = sum(old_pos.^2 .* wf.HOshape)
    #return exp.(wf.alpha * (old_r2 - r2))
    temp = old_pos.^2
    temp .= temp .- positions[:, p1].^2
    temp .= temp .* wf.HOshape
    return exp(wf.alpha * sum(temp))
end

function kinetic(positions, wf::SimpleGaussian, temp)::Float64
    temp .= positions.^2
    temp .= temp .* wf.HOshape.^2
    temp .= 2 * wf.alpha .* temp
    temp .= wf.HOshape .- temp
    return wf.alpha * sum(temp)
end

function QF(positions, p1, wf::SimpleGaussian)::Array{Float64}
    return -4 * wf.alpha .* positions[:, p1] .* wf.HOshape
end

function paramDer(positions, wf::SimpleGaussian, temp)::Float64
    temp .= positions.^2 .* wf.HOshape
    return -sum(temp) / 10
end