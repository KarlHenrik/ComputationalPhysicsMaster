struct Particles
    positions::Array{Float64}
    pos2::Array{Float64}
    pos_temp::Array{Float64}
    num::Int64
    dims::Int64
    Particles(positions) = new(positions, positions.^2, copy(positions), size(positions, 2), size(positions, 1))
end

function update_particles!(particles, p1, wf)
    
end