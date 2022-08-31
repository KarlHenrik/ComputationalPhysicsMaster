function consider!(walker::Walker{S, M}, wf::Slater, new_idx::Int64, move) where S where M <: Metro_Muts
    (; positions, metro_muts) = walker
    metro_muts.old_pos = positions[new_idx]    
    positions[new_idx] += move
    
    walker.new_amp = amplitude(wf, positions, new_idx)
    ratio = walker.new_amp / walker.old_amp
    
    return ratio
end

function accept!(walker, wf::Slater, new_idx, ham)
    (; samp_muts, positions) = walker
    walker.old_amp = wf_m.new_amp
    walker.accepted = true
    
    update_sample!(samp_muts, positions, wf, ham)
    return walker
end

# Slater needs to reset some columns
function deny!(walker, wf, new_idx)
    walker.positions[new_idx] .= walker.metro_m.old_pos
    walker.accepted = false
    return walker
end


"""
function consider!(walker::Walker{W, M, S}, system, new_idx) where M <: Metro_m_imp_slater where S <: Sample_m_energy_slater
    # Do the triple slater evaluation
    
    return walker
end

function consider!(walker::Walker{W, M, S}, system, new_idx) where M <: Metro_m_imp_slater
    # Do the single + QF slater evaluation
    
    return walker
end

function consider!(walker::Walker{W, M, S}, system, new_idx) where S <: Sample_m_energy_slater
    # Do the single + dder slater evaluation
    
    return walker
end
""";