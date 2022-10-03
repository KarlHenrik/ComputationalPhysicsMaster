abstract type Gaussian <: WaveFunction end

function paramDerHolder(wf::Gaussian)
    return 0.0
end

function paramDer!(samp_muts, positions, wf::Gaussian)
    pos_sum = 0.0
    for pos in positions
        pos_sum += pos^2
    end
    
    samp_muts.paramDer = -pos_sum / wf.n
    return samp_muts
end

function consider!(wf::Gaussian, positions, new_idx, old_pos)
    ratio = ratio_direct(wf, positions, new_idx, old_pos)
    return ratio
end

function consider_qf!(wf::Gaussian, positions, new_idx, old_pos)
    ratio = ratio_direct(wf, positions, new_idx, old_pos)
    newQF = QF(positions, new_idx, wf)
    
    return ratio, newQF
end

function accept!(wf::Gaussian, new_idx, new_pos)
    return wf
end