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

function accept!(walker, wf::NeuralNetwork, ham)
    (; samp_muts, position) = walker
    walker.old_amp = wf_m.new_amp
    walker.accepted = true
    
    update_sample!(samp_muts, positions, wf, ham)
    wf.QF_all_old .= wf.QF_all
    
    return walker
end