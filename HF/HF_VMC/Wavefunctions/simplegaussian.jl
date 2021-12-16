struct SimpleGaussian <: WaveFunction
    alpha::Float64
    HOshape::Vector{Float64}
    HOshape2::Vector{Float64}
    function SimpleGaussian(alpha, HOshape)
        return new(alpha, HOshape, HOshape.^2)
    end
end


function evaluate(particles, wf::SimpleGaussian)::Float64
    """
    Old wavefunc value term: exp(-alpha * old_r2)
    New wavefunc value term: exp(-alpha * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    exponent = 0.0
    for d in 1:particles.dims
        for i in d:particles.dims:length(particles.positions)
            exponent += particles.positions[i]^2 * wf.HOshape[d]
        end
    end
    return exp(-wf.alpha * exponent)
end

function kinetic(particles, wf::SimpleGaussian)::Float64
    kin_sum = 0.0
    for d in 1:particles.dims
        for i in d:particles.dims:length(particles.positions)
            kin_sum += particles.positions[i]^2 * wf.HOshape2[d]
        end
    end
    return wf.alpha * (sum(wf.HOshape) * particles.num - 2.0 * wf.alpha * kin_sum)
end

function QF(particles, idx, wf::SimpleGaussian) #maybe faster to put return type here?
    @views qf = -4.0 * wf.alpha .* particles.positions[idx] .* wf.HOshape
    return qf
end

function QF!(qf, particles, idx, wf::SimpleGaussian) #maybe faster to put return type here?
    @views qf .= -4.0 * wf.alpha .* particles.positions[idx] .* wf.HOshape
    return qf
end

function paramDer(particles, wf::SimpleGaussian)::Float64
    pos_sum = 0.0
    for d in 1:particles.dims
        for i in d:particles.dims:length(particles.positions)
            pos_sum += particles.positions[i]^2 * wf.HOshape[d]
        end
    end
    return -pos_sum / particles.num
end

function applyGradient(wf::SimpleGaussian, grad)
    return SimpleGaussian(wf.alpha - grad, wf.HOshape)
end