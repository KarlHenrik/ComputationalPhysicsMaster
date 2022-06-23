function consider!(walker, system, new_idx::Int64, move)
    start_idx = new_idx - (new_idx-1)%wf.dims
    dim = 
    
    new_idx = start_idx:start_idx + wf.dims - 1
    new_pos = walker.positions[new_idx]
    return ratio(walker.positions, new_idx, new_pos, wf)
end

function consider!(walker::Walker{W, M, S}, system, new_idx, move) where W where M <: Metro_m_imp where S
    walker.metro_m.old_pos .= walker.positions[new_idx]
    walker.metro_m.new_pos .= walker.positions[new_idx] .+ move
    
    ratio_ = ratio(walker.positions, new_idx, walker.metro_m.new_pos, wf)
    walker.positions[new_idx] .= walker.metro_m.new_pos
    QF!(walker.metro_m.newQF, walker.positions, new_idx, wf)
    return ratio_
end

function accept!(walker, system, new_idx, move)
    walker.positions[new_idx] += move
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