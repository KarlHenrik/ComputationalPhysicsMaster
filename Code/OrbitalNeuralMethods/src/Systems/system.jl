abstract type System end

include("basis.jl")

include("Integrals/interactions.jl")
include("Integrals/spatialintegrals.jl")
include("Integrals/spinintegrals.jl")

struct SpatialSystem{T} <: System
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    
    h::Array{Float64, 2} # one-body integrals
    u::Array{Float64, 4} # two-body integrals
    spfs::Vector{Vector{Float64}} # the basis functions evaluated on the grid
    
    grid::Vector{Float64}
    basis::T
    transform::Matrix{Float64}
    V::Interaction
end

function System(n, basis::SpatialBasis, grid, V::Interaction)
    l = basis.l
    spfs = spatial(basis, grid) # The basis functions evaluated on the grid
    h = onebody(basis, grid) # One body integrals
    u = twobody(basis, grid, V) # Two body integrals
    u .= u .- permutedims(u, [1, 2, 4, 3]) # Anti-symmetrizing u

    transform = la.I(l)
    return SpatialSystem{typeof(basis)}(n, l, h, u, spfs, grid, basis, transform, V)
end

function sp_energies(system)
    (; l, h) = system
    
    ϵ = zeros(l)
    @inbounds for q in 1:l
        ϵ[q] = h[q, q]
    end
    return ϵ
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

function sp_energies(system::PairingSystem)
    (; l, h, u) = system
    
    ϵ = zeros(l)
    @inbounds for q in 1:l
        ϵ[q] = h[q, q]
        for i in 1:n
            ϵ[q] += u[q, i, q, i] # I don't know why this is the setup in the book benchmark. See p.364 of An AdvancedCourse in Computational Nuclear Physics
        end
    end
    return ϵ
end

include("Integrals/transform.jl")
;