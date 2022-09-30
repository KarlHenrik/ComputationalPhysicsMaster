struct SlaterNN{N} <: WaveFunction
    n::Int64

    slater::Slater
    nn::NeuralNetwork{N}

    minisort::MiniSort
    nn_dder::Vector{Float64}
    nn_der::Vector{Float64}
    permtemp::Vector{Float64}
end
function SlaterNN(slater::Slater, nn::NeuralNetwork{N}) where N
    n = slater.n
    nn_dder, nn_der, permtemp = zeros(n), zeros(n), zeros(n)
    minisort = MiniSort(zeros(n))

    return SlaterNN{N}(n, slater, nn, minisort, nn_dder, nn_der, permtemp)
end

function private_wf(wf::SlaterNN{N}, positions) where N
    (; n, slater, nn) = wf
    slater = private_wf(slater, positions)
    minisort = MiniSort(positions)
    nn = private_wf(nn, minisort.x_sort)
    
    nn_dder, nn_der, permtemp = zeros(n), zeros(n), zeros(n)

    return SlaterNN{N}(n, slater, nn, minisort, nn_dder, nn_der, permtemp)
end


#* Do I keep two copies of sorted positions and reverse sorting indexes, or just one which I update?
#* I need the reverse sorting only for kinetic. I need the old sorting for not accepting.
#* Idea: Keep an old and new sorted. Only update reverse sorting when accepting. Yeah seems smart.

function consider_qf!(wf::SlaterNN, positions, new_idx::Int64, old_pos)
    (; slater, nn, minisort) = wf
    new_pos = positions[new_idx]
    sorted_pos, sort_idx = try_sort!(minisort, new_idx, new_pos)

    nn_ratio, nn_newQF = consider_qf!(nn, sorted_pos, sort_idx, old_pos)
    slater_ratio, slater_newQF = consider_qf!(slater, positions, new_idx, old_pos)
    
    newQF = slater_newQF + nn_newQF
    ratio = slater_ratio * nn_ratio

    return ratio, newQF
end

function consider!(wf::Slater, positions, new_idx::Int64, old_pos)
    (; slater, nn, minisort) = wf
    new_pos = positions[new_idx]
    sorted_pos, sort_idx = try_sort!(minisort, new_idx, new_pos)

    nn_ratio = consider!(nn, sorted_pos, sort_idx, old_pos)
    slater_ratio = consider!(slater, positions, new_idx, old_pos)
    
    ratio = slater_ratio * nn_ratio

    return ratio
end

function accept!(wf::SlaterNN, new_idx, new_pos)
    (; slater, nn, minisort) = wf
    accept!(slater, new_idx, new_pos)
    accept!(nn, new_idx, new_pos) # This is the wrong index for the nn, but it does not matter, as it's not used
    update_sort!(minisort, new_idx, new_pos) # Updating the sorted positions
    return wf
end

#* Only used for the old qf
function QF(positions, new_idx::Int64, wf::SlaterNN)
    (; slater, nn, minisort) = wf
    
    (; x_sort, p) = minisort

    return QF(positions, new_idx, slater) + QF(x_sort, p[new_idx], nn)
end

function permute_noalloc!(v, p, temp)
    temp .= v
    @inbounds for (i, idx) in enumerate(p)
        v[i] = temp[idx]
    end
    return v
end

function kinetic(positions, wf::SlaterNN)::Float64
    (; slater, nn, minisort, nn_dder, nn_der, permtemp) = wf
    (; x_sort, p_rev) = minisort
    
    nn_dder = dder!(nn_dder, x_sort, nn) # Double derivative with sorted positions
    permute_noalloc!(nn_dder, p_rev, permtemp) # Permuting double derivatives to unsored indexes
    
    nn_der .= 0.5 .* nn.QF_all_old # Derivative with sorted positions
    permute_noalloc!(nn_der, p_rev, permtemp) # Permuting derivatives to unsored indexes
    
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

function applyGradient(wf::SlaterNN, grad)
    nn = applyGradient(wf.nn, grad)

    wf = SlaterNN(wf.slater, nn)
    return wf
end

function paramDerHolder(wf::SlaterNN)
    return paramDerHolder(wf.nn)
end

function paramDer!(layer_grads::Vector{Dense_Grad}, positions, wf::SlaterNN)
    return paramDer!(layer_grads, positions, wf.nn)
end