function consider!(walker, system, new_idx::Int64, move)
    start_idx = new_idx - (new_idx-1)%wf.dims
    dim = 
    
    new_idx = start_idx:start_idx + wf.dims - 1
    new_pos = walker.positions[new_idx]
    return ratio(walker.positions, new_idx, new_pos, wf)
end

# 
function consider!(walker::Walker{W, M, S}, wf, new_idx, move) where W where M <: Metro_m_imp where S
    (; positions, metro_m, wf_m) = walker
    metro_m.old_pos .= positions[new_idx]
    
    positions[new_idx] .= positions[new_idx] .+ move
    
    wf_m.new_amp = amplitude(wf, positions)
    QF!(metro_m.newQF, wf, positions, new_idx)
    
    ratio = wf_m.new_amp / wf_m.old_amp
    return ratio, metro_m.newQF
end

function accept!(walker)
    walker.wf_m.old_amp = wf_m.new_amp
    walker.accepted = true
    return walker
end

function deny!(walker, new_idx)
    walker.positions[new_idx] .= walker.metro_m.old_pos
    walker.accepted = false
    return walker
end

function accept!(walker::Walker{W, M, S}, system) where W where M where S <: Sample_m_blocking
    #walker.positions[walker.metro_m.new_idx] = walker.metro_m.new_pos
    walker.sample_m.kinetic = kinetic(walker.positions, system.wf)
    return walker
end

function accept!(walker::Walker{W, M, S}, system) where W where M where S <: Sample_m_gradient
    #walker.positions[walker.metro_m.new_idx] = walker.metro_m.new_pos
    walker.sample_m.kinetic = kinetic(walker.positions, system.wf)
    walker.sample_m.paramDer = paramDer(walker.positions, system.wf)
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