struct SlaterNN <: WaveFunction
    n::Int64

    slater::Slater
    nn::NeuralNetwork
end

function private_wf(wf::SlaterNN, positions)
    (; slater, nn) = wf
    slater = private_wf(slater, positions)
    nn = private_wf(nn, positions)

    return SlaterNN(n, slater, nn)
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
    (; slater, nn) = wf
    p = sortperm(positions)
    return QF(positions, new_idx, slater) + QF(positions[p], new_idx[p], nn)
end

function kinetic(positions, wf::SlaterNN)::Float64
    (; slater, nn) = wf
    p = sortperm(positions)
    p_back = sortperm(p)
    nn_dder = dder(positions[p], nn)
    nn_dder .= nn_dder[p_back]
    nn_der = 0.5 .* nn.QF_all_old
    nn_der .= nn_der[p_back]

    (; amp_up, amp_down, kin_mat, n) = slater
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