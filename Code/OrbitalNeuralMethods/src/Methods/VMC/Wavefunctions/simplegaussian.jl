struct SimpleGaussian <: WaveFunction
    α::Float64
    dims::Int64
    num::Int64
    HOshape::Vector{Float64}
    HOshape2::Vector{Float64}
    function SimpleGaussian(dims, num; α, HOshape=ones(dims))
        @assert dims == length(HOshape)
        return new(α, dims, num, HOshape, HOshape.^2)
    end
end
private_wf(wf::SimpleGaussian) = SimpleGaussian(wf.dims, wf.num, α=wf.α, HOshape=wf.HOshape)

# Index from metro does not cover all dimentions
function ratio_direct(wf::SimpleGaussian, positions, new_idx::Int64, old_pos)
    """
    Old wavefunc value term: exp(-α * old_r2)
    New wavefunc value term: exp(-α * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    d = (new_idx-1)%wf.dims +1
    ratio_sum = (old_pos^2 - positions[new_idx]^2) * wf.HOshape2[d]
    return exp(wf.α * ratio_sum)
end

# Index from importance covers all dims
function ratio_direct(wf::SimpleGaussian, positions, new_idx, old_pos)
    """
    Old wavefunc value term: exp(-α * old_r2)
    New wavefunc value term: exp(-α * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    ratio_sum = 0.0
    for d in 1:wf.dims
        ratio_sum += (old_pos[d]^2 - positions[new_idx[d]]^2) * wf.HOshape2[d]
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

function QF!(qf, positions, idx, wf::SimpleGaussian)
    @views qf .= -4.0 * wf.α .* positions[idx] .* wf.HOshape
    return qf
end

function paramDer!(samp_muts, positions, wf::SimpleGaussian)
    pos_sum = 0.0
    for d in 1:wf.dims
        for i in d:wf.dims:length(positions)
            pos_sum += positions[i]^2 * wf.HOshape[d]
        end
    end
    samp_muts.paramDer = -pos_sum / wf.num
    return samp_muts
end

function paramDerHolder(wf::SimpleGaussian)
    return 0.0
end

function applyGradient(wf::SimpleGaussian, grad)
    return SimpleGaussian(wf.dims, wf.num, α = wf.α - grad, HOshape=wf.HOshape)
end