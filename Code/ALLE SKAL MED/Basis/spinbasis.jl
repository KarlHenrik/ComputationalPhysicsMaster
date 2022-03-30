struct SpinBasis{T} <: SpatialBasis
    base::T # The basis with no spin
    l::Int64
    function SpinBasis(base::SpatialBasis)
        return new{typeof(base)}(base, 2 * base.l)
    end
end

function spatial(basis::SpinBasis, grid)
    n = length(grid)
    l = basis.l
    res = [zeros(n) for i in 1:l]
    
    nospin = zeros(l÷2)
    
    @inbounds for i in 1:n
        evaluate!(nospin, grid[i], basis.base) # the basis functions evaluated at x
        for j in 1:l÷2
            res[2j-1][i] = nospin[j]
            res[2j][i] = nospin[j]
        end
    end
    
    return res
end

function onebody(basis::SpinBasis, grid)
    h = onebody(basis.base, grid)
    h = kron(h, [1 0; 0 1])
    return h
end

function twobody(basis::SpinBasis, grid, a)
    u = twobody(basis.base, grid, a)
    u = add_spin_u(u)
    return u
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