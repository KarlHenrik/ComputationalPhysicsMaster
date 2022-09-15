struct SimpleGaussian <: Gaussian
    α::Float64
    n::Int64

    function SimpleGaussian(n; α)
        return new(α, n)
    end
end
private_wf(wf::SimpleGaussian, positions) = SimpleGaussian(wf.n, α=wf.α)

# Index from metro does not cover all dimentions
function ratio_direct(wf::SimpleGaussian, positions, new_idx::Int64, old_pos)
    """
    Old wavefunc value term: exp(-α * old_r2)
    New wavefunc value term: exp(-α * r2)
    All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
    """
    ratio_sum = (old_pos^2 - positions[new_idx]^2)
    return exp(wf.α * ratio_sum)
end

function kinetic(positions, wf::SimpleGaussian)::Float64
    kin_sum = 0.0
    for pos in positions
        kin_sum += pos^2
    end

    return wf.α * (wf.n - 2.0 * wf.α * kin_sum)
end

function QF(positions, idx, wf::SimpleGaussian)
    qf = -4.0 * wf.α .* positions[idx]
    return qf
end