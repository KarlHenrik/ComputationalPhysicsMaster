struct SlaterNN{T} <: WaveFunction
    n::Int64

    slater::Slater
    nn::NeuralNetwork{T}

    p::Vector{Int64}
    p_back::Vector{Int64}
    nn_dder::Vector{Float64}
    nn_der::Vector{Float64}
    sorttemp::Vector{Float64}
end
function SlaterNN(slater::Slater, nn::NeuralNetwork)
    n = slater.n
    p, p_back, nn_dder, nn_der, sorttemp = zeros(Int, n), zeros(Int, n), zeros(n), zeros(n), zeros(n)

    return SlaterNN{typeof(nn.hes_jac_config)}(n, slater, nn, p, p_back, nn_dder, nn_der, sorttemp)
end

function private_wf(wf::SlaterNN, positions)
    (; slater, nn) = wf
    slater = private_wf(slater, positions)
    nn = private_wf(nn, positions)

    return SlaterNN(slater, nn)
end

function consider_qf!(wf::SlaterNN, positions, new_idx::Int64, old_pos)
    (; slater, nn) = wf
    
    p = sortperm(positions)
    nn_ratio, nn_newQF = consider_qf!(nn, positions[p], p[new_idx], old_pos)
    slater_ratio, slater_newQF = consider_qf!(slater, positions, new_idx, old_pos)
    
    newQF = slater_newQF + nn_newQF
    ratio = slater_ratio * nn_ratio

    return ratio, newQF
end

function accept!(wf::SlaterNN, new_idx)
    (; slater, nn) = wf
    accept!(slater, new_idx)
    accept!(nn, new_idx) # This is the wrong index for the nn, but it does not matter, as it's not used
    return wf
end

#* Only used for the old qf
function QF(positions, new_idx::Int64, wf::SlaterNN)
    (; slater, nn, sorttemp) = wf
    p = sortperm(positions)

    @inbounds for (i, idx) in enumerate(p)
        sorttemp[i] = positions[idx]
    end
    
    return QF(positions, new_idx, slater) + QF(sorttemp, p[new_idx], nn)
end

function permute_noalloc!(v, p, temp)
    temp .= v
    @inbounds for (i, idx) in enumerate(p)
        v[i] = temp[idx]
    end
    return v
end

function kinetic(positions, wf::SlaterNN)::Float64
    (; slater, nn, p, p_back, nn_dder, nn_der, sorttemp) = wf
    p = sortperm!(p, positions)
    p_back = sortperm!(p_back, p)
    
    for (i, idx) in enumerate(p)
        sorttemp[i] = positions[idx]
    end
    nn_dder = dder!(nn_dder, sorttemp, nn)

    #@inbounds nn_dder .= nn_dder[p_back]
    permute_noalloc!(nn_dder, p_back, sorttemp)
    
    nn_der .= 0.5 .* nn.QF_all_old
    #@inbounds nn_der .= nn_der[p_back]
    permute_noalloc!(nn_der, p_back, sorttemp)
    
    (; amp_up, amp_down, der_mat, kin_mat, n) = slater
    kin = 0.0
    #TODO maybe find faster way of doing this than slices
    for particle in 1:n÷2
        @views slater_dder = ratio_new_old_det(amp_up, kin_mat[particle, :], particle)
        @views slater_der = ratio_new_old_det(amp_up, der_mat[particle, :], particle)
        kin += slater_dder + nn_dder[particle] + 2 * slater_der * nn_der[particle]
    end
    for (i, particle) in enumerate(n÷2+1:n)
        @views slater_dder = ratio_new_old_det(amp_down, kin_mat[particle, :], i)
        @views slater_der = ratio_new_old_det(amp_down, der_mat[particle, :], i)
        kin += slater_dder + nn_dder[particle] + 2 * slater_der * nn_der[particle]
    end

    return -0.5 * kin
end