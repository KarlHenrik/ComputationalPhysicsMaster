function consider!(state, system, new_idx::Int64, move)
    start_idx = new_idx - (new_idx-1)%wf.dims
    dim = 
    
    new_idx = start_idx:start_idx + wf.dims - 1
    new_pos = state.positions[new_idx]
    return ratio(state.positions, new_idx, new_pos, wf)
end

function consider!(state::State{W, M, S}, system, new_idx, move) where W where M <: Metro_m_imp where S
    state.metro_m.old_pos .= state.positions[new_idx]
    state.metro_m.new_pos .= state.positions[new_idx] .+ move
    
    ratio_ = ratio(state.positions, new_idx, state.metro_m.new_pos, wf)
    state.positions[new_idx] .= state.metro_m.new_pos
    QF!(state.metro_m.newQF, state.positions, new_idx, wf)
    return ratio_
end

function accept!(state, system, new_idx, move)
    state.positions[new_idx] += move
    return state
end

function accept!(state::State{W, M, S}, system) where W where M where S <: Sample_m_blocking
    #state.positions[state.metro_m.new_idx] = state.metro_m.new_pos
    state.sample_m.kinetic = kinetic(state.positions, system.wf)
    return state
end

function accept!(state::State{W, M, S}, system) where W where M where S <: Sample_m_gradient
    #state.positions[state.metro_m.new_idx] = state.metro_m.new_pos
    state.sample_m.kinetic = kinetic(state.positions, system.wf)
    state.sample_m.paramDer = paramDer(state.positions, system.wf)
    return state
end


"""
function consider!(state::State{W, M, S}, system, new_idx) where M <: Metro_m_imp_slater where S <: Sample_m_energy_slater
    # Do the triple slater evaluation
    
    return state
end

function consider!(state::State{W, M, S}, system, new_idx) where M <: Metro_m_imp_slater
    # Do the single + QF slater evaluation
    
    return state
end

function consider!(state::State{W, M, S}, system, new_idx) where S <: Sample_m_energy_slater
    # Do the single + dder slater evaluation
    
    return state
end
""";