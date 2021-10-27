struct Particles
    positions::Vector{Float64}
    dims::Int64
    num::Int64
    
    # Temporary vectors used for importance sampling
    temp_vec::Vector{Float64}
    temp_vec2::Vector{Float64}
    temp_vec3::Vector{Float64}
end

function Particles(dims, num, rng)
    positions = zeros(dims * num)
    for i in 1:dims * num
        positions[i] = (Random.rand(rng, Float64) .- 0.5) .* 2
    end
    temp_vec = zeros(dims)
    temp_vec2 = zeros(dims)
    temp_vec3 = zeros(dims)
    return Particles(positions, dims, num, temp_vec, temp_vec2, temp_vec3)
end

function Particles(dims, num, rng, wf::WaveFunction)
    return Particles(dims, num, rng)
end

#=
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
            
            candidate = (Random.rand(rng, Float64, dims) .- 0.5) .* 2
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
    positions = sa.SizedVector{dims * num}(reduce(vcat, positions))
    return Particles(positions, dims, num, temp_vec)
end
=#