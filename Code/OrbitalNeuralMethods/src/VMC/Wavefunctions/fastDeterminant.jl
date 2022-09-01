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