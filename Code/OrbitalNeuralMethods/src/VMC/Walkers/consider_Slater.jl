# ------------------- No first or second derivative ------------------------
function consider!(walker::Walker{S, M}, wf::Slater, new_idx::Int64, move::Float64) where S where M <: Metro_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos = positions[new_idx]    
    positions[new_idx] += move
    
    computeNewRows!(wf, walker, positions[new_idx], new_idx)

    ratio = ratio_direct(wf)

    return ratio
end

function consider!(newQF, walker::Walker{S, M}, wf::Slater, new_idx::Int64, move::Float64) where S where M <: Imp_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos = positions[new_idx]    
    positions[new_idx] += move
    
    wf = computeNewRows!(wf, walker, positions[new_idx], new_idx)
    (; amp_up, amp_down, new_der) = wf

    ratio = ratio_direct(wf)
    if new_idx < num÷2
        newQF .+= ratio_new_old_det(amp_up, new_der, idx)
    else
        newQF .+= ratio_new_old_det(amp_down, new_der, idx-num÷2)
    end

    return ratio, newQF
end

function accept!(walker::Walker{S, M}, wf::Slater, new_idx, ham) where S <: No_Muts where M <: Metro_Muts
    (; samp_muts) = walker
    walker.accepted = true
    
    wf = setNewRows!(wf, walker, new_idx)

    update_sample!(samp_muts, positions, wf, ham)
    return walker
end


;