function update_sample!(samp_muts::No_Muts, positions, wf, ham)
    return samp_muts
end

function update_sample!(samp_muts::E_Muts, positions, wf, ham)
    samp_muts.kinetic = kinetic(positions, wf)
    samp_muts.potential = potential(positions, ham)
    samp_muts.E = samp_muts.kinetic + samp_muts.potential
    samp_muts.E2 = samp_muts.E^2
    
    return samp_muts
end

function update_sample!(samp_muts::Grad_Muts, positions, wf, ham)
    samp_muts.kinetic = kinetic(positions, wf)
    samp_muts.potential = potential(positions, ham)
    samp_muts.E = samp_muts.kinetic + samp_muts.potential
    samp_muts.E2 = samp_muts.E^2
    
    samp_muts = paramDer!(samp_muts, positions, wf)
    
    return samp_muts
end


# ------------------ Denying moves -------------------

# Importance + Anything but Slater
function deny!(walker, wf, new_idx::UnitRange{Int64})
    walker.positions[new_idx] .= walker.metro_muts.old_pos
    walker.accepted = false
    return walker
end

# Metropolis + Anything but Slater
function deny!(walker, wf, new_idx::Int64)
    walker.positions[new_idx] = walker.metro_muts.old_pos
    walker.accepted = false
    return walker
end