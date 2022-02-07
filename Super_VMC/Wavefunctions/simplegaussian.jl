struct SimpleGaussian <: WaveFunction
    α::Float64
    dims::Int64
    num::Int64
    HOshape::Vector{Float64}
    HOshape2::Vector{Float64}
    function SimpleGaussian(HOshape, num, α)
        return new(α, length(HOshape), num, HOshape, HOshape.^2)
    end
end
SimpleGaussian(HOshape, num; α) = SimpleGaussian(HOshape, num, α)
SimpleGaussian(dims::Int64, num; α) = SimpleGaussian(ones(dims), num, α)

function ratio(dist, positions, new_idx, new_pos, wf::SimpleGaussian)::Float64
    """
    Old wavefunc value term: exp(-α * old_r2)
    New wavefunc value term: exp(-α * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    start_idx = new_idx - (new_idx-1)%wf.dims - 1
    ratio_sum = 0.0
    for d in 1:wf.dims
        ratio_sum += (positions[start_idx+d]^2 - new_pos[d]^2) * wf.HOshape2[d]
    end
    return exp(wf.α * ratio_sum)
end


function kinetic(positions, wf::SimpleGaussian)::Float64
    kin_sum = 0.0
    for d in 1:wf.dims
        for i in d:wf.dims:length(positions)
            kin_sum += positions[i]^2 * wf.HOshape2[d]
        end
    end
    return wf.α * (sum(wf.HOshape) * wf.num - 2.0 * wf.α * kin_sum)
end

function QF(positions, idx, wf::SimpleGaussian)
    @views qf = -4.0 * wf.α .* positions[idx] .* wf.HOshape
    return qf
end

function QF!(qf, positions, idx, wf::SimpleGaussian)
    @views qf .= -4.0 * wf.α .* positions[idx] .* wf.HOshape
    return qf
end

function paramDer(positions, wf::SimpleGaussian)::Float64
    pos_sum = 0.0
    for d in 1:wf.dims
        for i in d:wf.dims:length(positions)
            pos_sum += positions[i]^2 * wf.HOshape[d]
        end
    end
    return -pos_sum / wf.num
end

function applyGradient(wf::SimpleGaussian, grad)
    return SimpleGaussian(wf.HOshape, num, wf.α - grad)
end