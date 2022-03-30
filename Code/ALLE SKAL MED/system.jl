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
    # The basis functions evaluated on the grid
    spfs = spatial(basis, grid)

    # One body integrals
    h = onebody(basis, grid)

    # Two body integrals
    inner = inner_ints(spfs, grid, a)
    u = outer_int(spfs, grid, inner)

    # Adding spin
    l = 2 * basis.l
    h = kron(h, [1 0; 0 1])
    u = add_spin_u(u)
    spfs = [spfs[(i + 1)÷2] for i in 1:l]

    # Anti-symmetrizing u
    u .= u .- permutedims(u, [1, 2, 4, 3])

    return SpatialSystem{typeof(basis)}(n, l, a, h, u, spfs, grid, basis)
end

function trapz(f_vals, grid)
    val = sum(f_vals)
    val = val - 0.5 * (f_vals[1] + f_vals[end])

    return val * (grid[2] - grid[1])
end

function inner_ints(spfs, grid, a)
    l = length(spfs)
    inner_int = zeros(l, l, length(spfs[1]))
    cs = [zero(grid) for i in 1:Threads.nthreads()]
    fs = [zero(grid) for i in 1:Threads.nthreads()]
    
    @inbounds Threads.@threads for xi in eachindex(grid)
        x1 = grid[xi]
        f_vals = fs[Threads.threadid()]
        coulomb = cs[Threads.threadid()]
        
        coulomb .= 1 ./ sqrt.( (grid .- x1).^2 .+ a.^2 )
        for κ in 1:l
            for λ in κ:l
                f_vals .= conj.(spfs[κ]) .* coulomb .* spfs[λ]
                res = trapz(f_vals, grid)
                inner_int[κ, λ, xi] = res
                inner_int[λ, κ, xi] = res
            end
        end
    end
    return inner_int
end


function outer_int(spfs, grid, inner_ints)
    l = length(spfs)
    outer_int = zeros(l, l, l, l)
    
    fs = [zero(grid) for i in 1:Threads.nthreads()]
    is = [zero(grid) for i in 1:Threads.nthreads()]

    @inbounds Threads.@threads for κ in 1:l
        f_vals = fs[Threads.threadid()]
        inner = is[Threads.threadid()]
        
        for λ in 1:l
            @views inner .= inner_ints[κ, λ, :]
            for μ in 1:l
                for ν in 1:l
                    f_vals .= conj.(spfs[μ]) .* inner .* spfs[ν]
                    outer_int[μ, κ, ν, λ] = trapz(f_vals, grid)
                end
            end
        end
    end
    return outer_int
end

function add_spin_u(u_old)
    """
    When we make a spin-up and spin-down duplicate of each basis function, we get many new integrals to compute
    between these new basis functions. However, most of these integrals will be zero due to opposite spins between
    the basis functions, or they will be the same as the old integrals, if the spins align.
    
    We loop over the latter two indices in our two-body integral, ν and λ. The matrix of elements at the indeces
    u_new[:, :, ν, λ] then correspond to integrals where the first two basis functions in the integral have spins
    alternating up and down.
    
    'up-up'   'up-down'   up-up   ...
    'down-up' 'down-down' down-up ...
     up-up     up-down    up-up   ...
     ...       ...        ...     ...
    
    This 2x2 pattern of spins (that we have marked with '') repeats, and all elements in the 2x2 pattern has the same
    set of spatial functions. But only one element will be non-zero, the one that has the same spins as the basis
    functions numbered ν and λ.
    
    The code below finds which one of these elements in the 2x2 pattern will be non-zero, and then uses this to turn
    the two-body-integrals without spin into the two-body-integrals with spin.
    """
    l = size(u_old)[1] * 2
    u_new = zeros((l, l, l, l))
    
    for ν in 1:l
        for λ in 1:l
            if (ν % 2 == 1)     # --- UP ---
                if (λ % 2 == 1) # UP - UP
                    pattern = [1 0; 0 0]
                else            # UP - DOWN
                    pattern = [0 1; 0 0]
                end
            else                # --- DOWN ---
                if (λ % 2 == 1) # DOWN - UP
                    pattern = [0 0; 1 0]
                else            # DOWN - DOWN
                    pattern = [0 0; 0 1]
                end
            end
            
            ν_old = Int( ceil(ν / 2) )
            λ_old = Int( ceil(λ / 2) )
            @views kron!(u_new[:, :, ν, λ], u_old[:, :, ν_old, λ_old], pattern)
        end
    end
    return u_new
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