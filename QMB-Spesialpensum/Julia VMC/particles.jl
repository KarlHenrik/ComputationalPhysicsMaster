import Random
import LinearAlgebra as la

include("Wavefunctions/wavefunctions.jl")

struct Particles{M, V}
    positions::M
    dims::Int64
    num::Int64
    temp_vec::V
end

function Particles(dims, num, rng, wf::WaveFunction)
    positions = Vector{sa.MVector{dims, Float64}}(undef, num)
    for new_particle in 1:num
        positions[new_particle] = sa.MVector{dims}( (rand(rng, Float64, dims) .- 0.5) .* 2 )
    end
    temp_vec = copy(positions[1])
    return Particles(positions, dims, num, temp_vec)
end

function Particles(dims, num, rng, wf::Correlated)
    positions = Vector{sa.MVector{dims, Float64}}(undef, num)
    for new_particle in 1:num
        notFound = true
        attempts = 0
        local candidate
        while notFound
            attempts += 1
            if attempts > 1000
                error("No more room for particles. Reduce the number of particles or the hard-shell radius")
            end
            
            candidate = sa.MVector{dims}( (rand(rng, Float64, dims) .- 0.5) .* 2 )
            notFound = false
            for placed_particle in 1:new_particle - 1
                if la.norm(positions[placed_particle] .- candidate) < wf.a
                    notFound = true
                end
            end
        end
        positions[new_particle] = candidate
    end
    
    temp_vec = copy(positions[1])
    return Particles(positions, dims, num, temp_vec)
end