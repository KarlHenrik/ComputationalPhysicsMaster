abstract type System end

struct SpatialSystem{T} <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    a::Float64 # Coulomb shielding
    
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
    spfs::Vector{Vector{Float64}} # the basis functions evaluated on the grid
    
    grid::Vector{Float64}
    basis::T
end

function System(n, basis::SpatialBasis, grid, a)
    l = basis.l
    spfs = spatial(basis, grid) # The basis functions evaluated on the grid
    h = onebody(basis, grid) # One body integrals
    u = twobody(basis, grid, a) # Two body integrals
    u .= u .- permutedims(u, [1, 2, 4, 3]) # Anti-symmetrizing u

    return SpatialSystem{typeof(basis)}(n, l, a, h, u, spfs, grid, basis)
end

struct PairingSystem <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
end

function System(n, basis::Pairing)
    (; l, states) = basis
    
    h = zeros((l, l))
    u = zeros((l, l, l, l))
    for p in 1:l
        for q in 1:l
            # One body integrals
            h[p, q] = Ĥ₀(states[p], states[q])

            for r in 1:l
                for s in 1:l
                    # Two body integrals
                    u[p, q, r, s] = V̂(states[p], states[q], states[r], states[s])
                end
            end
        end
    end

    return PairingSystem(n, l, h, u)
end
;