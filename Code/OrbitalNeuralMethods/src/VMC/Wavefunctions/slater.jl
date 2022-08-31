struct Slater{T} <: WaveFunction
    C::Matrix{Float64}
    l::Int64
    basis::T

    dims::Int64
    num::Int64

    new_row::Vector{Float64}

    amp_mat::Matrix{Float64}
    der_mat::Matrix{Float64}
    kin_mat::Matrix{Float64}

    function Slater(num, dims, C, basis)
        l = basis.l

        amp_mat = zero(C)
        der_mat = zero(C)
        kin_mat = zero(C)

        return new(C, l, basis, dims, num, amp_mat, der_mat, kin_mat)
    end
end

function Slater(state::HFState)
    (; C, system) = state
    (; basis, transform) = system
    C = transform * C
    return Slater(C, basis)
end

function Slater(system::System)
    (; basis, transform) = system
    C = transform
    return Slater(C, basis)
end

function evaluate!(wf::Slater{Basis}, x)
    (; basis_eval, basis) = wf
    basis_eval .= evaluate!(basis_eval, x, basis) # the basis functions evaluated at x
end

function evaluate!(wf::Slater{SpinBasis}, x)
    (; basis_eval, nospin, basis) = wf
    
    nospin .= evaluate!(nospin, x, basis.base) # the basis functions evaluated at x
    @inbounds for i in eachindex(nospin) # doubling the basis functions to include spin
        basis_eval[2i-1] = nospin[i]
        basis_eval[2i] = nospin[i]
    end
end

function amplitude(wf, positions, new_idx)
    (; C, l, basis_eval, amp_mat, num, new_row) = wf
    x = positions[new_idx]
    basis_eval = evaluate!(wf, x) # the basis functions evaluated at x
    
    for row in 1:num
        ϕ_col_x = 0.0
        for j in 1:l
            ϕ_col_x += C[j, row] * basis_eval[j]
        end
        new_row[i, row] = ϕ_col_x
    end
    
    return det_new_row(amp_mat, new_row, new_idx)
end

function kinetic(positions, wf::Slater)::Float64
    (; amp_mat, kin_mat, num) = wf
    kin = 0.0
    for row in 1:num
        kin += ratio_new_old_det(amp_mat, kin_mat[row, :], row)
    end
    return kin
end

# ------------------ Functinality for computing the change in determinant when changing only one row ---------------------

mutable struct Fast_Det
    const D::Matrix{Float64}
    const D_inv::Matrix{Float64}
    const n::Int64
    det::Float64

    
    function Fast_Det(D)
        D_inv = la.inv(D)
        det = la.det(D)
        n = size(D)[1]
        return new(copy(D), D_inv, n, det)
    end
end

function det(fd::Fast_Det)
    return fd.det
end

function det_new_row(fd::Fast_Det, new_row::Vector{Float64}, i)
    (; D_inv, n, det) = fd
    R = 0.0
    for j in 1:n
        R += new_row[j] * D_inv[j, i]
    end
    return det * R
end

function ratio_new_old_det(fd::Fast_Det, new_row::Vector{Float64}, i)
    (; D_inv, n) = fd
    R = 0.0
    for j in 1:n
        R += new_row[j] * D_inv[j, i]
    end
    return R
end

function change_row!(fd, new_row, i)
    (; D, D_inv, n) = fd
    R = 0.0
    for j in 1:n
        R += new_row[j] * D_inv[j, i]
    end
    fd.det = fd.det * R

    for j in 1:n
        if j != i
            S = 0.0
            for l in 1:n
                S += new_row[l] * D_inv[l, j]
            end
        
            for k in 1:n
                D_inv[k, j] = D_inv[k, j] - S / R * D_inv[k, i]
            end
        end
    end

    for k in 1:n
        D[i, k] = new_row[k]
        D_inv[k, i] /= R
    end
    
    return fd
end