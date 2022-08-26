# Metropolis + Gaussian/Correlated
function consider!(walker::Walker{S, M}, wf::Union{SimpleGaussian, Correlated}, new_idx, move) where S where M <: Metro_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos = positions[new_idx]
    
    positions[new_idx] = positions[new_idx] + move
    
    ratio = ratio_direct(wf, positions, new_idx, metro_muts.old_pos)
    return ratio
end

# Importance + Gaussian/Correlated
function consider!(newQF, walker::Walker{S, M}, wf::Union{SimpleGaussian, Correlated}, new_idx, move) where S where M <: Imp_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos .= positions[new_idx]
    
    positions[new_idx] .= positions[new_idx] .+ move
    
    ratio = ratio_direct(wf, positions, new_idx, metro_muts.old_pos)
    newQF = QF!(newQF, positions, new_idx, wf)
    
    return ratio, newQF
end

function accept!(walker, wf::Union{SimpleGaussian, Correlated}, new_idx, ham)
    (; samp_muts, positions) = walker
    walker.accepted = true
    
    update_sample!(samp_muts, positions, wf, ham)
    return walker
end