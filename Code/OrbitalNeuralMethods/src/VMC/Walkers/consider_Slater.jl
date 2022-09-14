# ------------------- No first or second derivative ------------------------
function consider!(walker::Walker{S, M}, wf::Slater, new_idx::Int64, move::Float64) where S where M <: Metro_Muts
    (; positions, metro_muts) = walker
    walker.old_pos = positions[new_idx]
    positions[new_idx] += move
    
    wf = computeNewRows!(wf, positions[new_idx])
    (; amp_up, amp_down, new_amp, n) = wf
    
    if new_idx <= n÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-n÷2)
    end

    return ratio
end

#TODO maybe look at making new_idx and move into ints and floats
function consider!(walker::Walker{S, M}, wf::Slater, new_idx::Int64, move::Float64) where S where M <: Imp_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos = positions[new_idx]
    positions[new_idx] += move
    
    wf = computeNewRows!(wf, positions[new_idx])
    (; amp_up, amp_down, new_amp, new_der, n) = wf
    
    if new_idx <= n÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
        newQF = ratio_new_old_det(amp_up, new_der, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-n÷2)
        newQF = ratio_new_old_det(amp_down, new_der, new_idx-n÷2)
    end

    return ratio, newQF
end

function accept!(walker::Walker{S, M}, wf::Slater, new_idx::Int64, ham) where S where M
    (; samp_muts, positions) = walker
    walker.accepted = true
    
    wf = setNewRows!(wf, new_idx)

    update_sample!(samp_muts, positions, wf, ham)
    return walker
end


;