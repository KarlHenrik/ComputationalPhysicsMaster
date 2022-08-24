# --------------- Considering moves ----------------------

# Metropolis + Gaussian/Correlated/NeuralNetwork + ParamDer is set up for NN
function consider!(walker::Walker{S, Q}, wf, new_idx, move) where S where Q <: No_muts
    (; positions, qf_muts) = walker
    qf_muts.old_pos .= positions[new_idx]
    
    positions[new_idx] .= positions[new_idx] .+ move
    
    walker.new_amp = amplitude(wf, positions)
    
    ratio = walker.new_amp / walker.old_amp
    return ratio
end

# Importance + Gaussian/Correlated
function consider!(newQF, walker::Walker{S, Q}, wf, new_idx, move) where S where Q <: QF_muts
    (; positions, qf_muts) = walker
    qf_muts.old_pos .= positions[new_idx]
    
    positions[new_idx] .= positions[new_idx] .+ move
    
    walker.new_amp = amplitude(wf, positions)
    newQF = QF!(newQF, wf, positions, new_idx)
    
    ratio = walker.new_amp / walker.old_amp
    return ratio, newQF
end

# Importance + NeuralNetwork + ParamDer is set up
function consider!(newQF, walker::Walker{S, Q}, wf::NeuralNetwork, new_idx, move) where S where Q <: QF_muts
    (; positions, qf_muts) = walker
    qf_muts.old_pos .= positions[new_idx]
    
    positions[new_idx] .= positions[new_idx] .+ move
    
    walker.new_amp = amplitude(wf, positions)
    qf_all!(wf)
    newQF .= wf.QF_all[new_idx]
    
    ratio = walker.new_amp / walker.old_amp
    return ratio, newQF
end

# --------------- Accepting moves --------------------

function accept!(walker, wf, ham)
    (; samp_muts, position) = walker
    walker.old_amp = wf_m.new_amp
    walker.accepted = true
    
    update_sample!(samp_muts, positions, wf, ham)
    update_wf_olds!(wf) #NN QF_all
    
    return walker
end

function update_sample!(samp_muts::E_muts, positions, wf, ham)
    samp_muts.kinetic = kinetic(positions, wf)
    samp_muts.potential = potential(positions, ham)
    samp_muts.E = samp_muts.kinetic + samp_muts.potential
    samp_muts.E2 = samp_muts.E^2
    
    return samp_muts
end

function update_sample!(samp_muts::Grad_muts, positions, wf, ham)
    samp_muts.kinetic = kinetic(positions, wf)
    samp_muts.potential = potential(positions, ham)
    samp_muts.E = samp_muts.kinetic + samp_muts.potential
    samp_muts.E2 = samp_muts.E^2
    
    samp_muts.paramDer = paramDer!(samp_muts.paramDer, wf)
    
    return samp_muts
end


# ------------------ Denying moves -------------------

function deny!(walker, new_idx)
    walker.positions[new_idx] .= walker.metro_m.old_pos
    walker.accepted = false
    return walker
end

# Slater needs to reset some columns
function deny!(walker, new_idx)
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