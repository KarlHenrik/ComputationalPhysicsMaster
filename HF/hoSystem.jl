struct HOFunc
    order::Float64
    factors::Vector{Float64}
    parity::Float64
    ω::Float64

    function HOFunc(n, ω)
        if n > 50
            print("HO functions past n=50 are not accurate enough. Use a recursive definition.")
        end
        factors = Float64[]
        
        ξfac = √ω
        
        for m in n÷2:-1:0
            factor = (-1)^m * (ω / π)^0.25 * ξfac^(n - 2m) * ho_factor(m, n)
            
            push!(factors, factor)
        end
        return new(n, factors, n % 2, ω)
    end
end

function ho_factor(m, n)
    # The harmonic oscillator basis functions need some factors of medium or very small size.
    # However, computing these factors requires factorials and exponentials that quickly 
    # get too large. By taking the logarithm of these factors, we don't need to
    # compute any large numbers, but can use logarithm tricks to make the computation easy.
    
    n_fac = sum([log(i) for i in 1:n])
    m_fac = sum([log(i) for i in 1:m])
    n_minus_2m_fac = sum([log(i) for i in 1:(n - 2m)])
    
              #√(      n!  /    2^(n)  ) *    2^(n - 2m)     / ( m!   *  (n - 2m)!    )
    log_fac = 0.5 * (n_fac - n * log(2)) + (n - 2m) * log(2) - (m_fac + n_minus_2m_fac)
    return exp(log_fac)
end

function compute(x, ho::HOFunc)
    result = zero(x)
    x2 = x^2
    xm = (x * ho.parity) + (1.0 - ho.parity) # xm = 1 if parity == 1 else xm = x
    for factor in ho.factors
        result += factor * xm
        xm *= x2
    end
    result = result * exp(-ho.ω * x^2 / 2)
    return result
end

function potential(x, ho::HOFunc)
    return 0.5 * ho.ω^2 * x^2
end

function E_n(ho::HOFunc)
    return (ho.order + 0.5) * ho.ω
end

# --------------------------------- The basis functions collected --------------------------------

abstract type Basis end

struct HOBasis <: Basis
    funcs::Vector{HOFunc}
    a::Float64
end

function HOBasis(l, ω, a)
    funcs = [HOFunc(i, ω) for i in 0:l-1]
    return HOBasis(funcs, a)
end

# --------------------------------- The system ----------------------------------

struct System{T<:Basis}
    n::Int64 # number of particles
    l::Int64   # number of basis functions
    
    h::Array{Float64, 2} # one-body integral
    u::Array{Float64, 4} # two-body integral
    spfs::Vector{Vector{Float64}} # the basis functions evaluated on the grid
    
    basis::T
    grid::Vector{Float64}
    
    function System(n, basis, grid)
        # The basis functions evaluated on the grid
        spfs = [compute.(grid, (ψ_i,)) for ψ_i in basis.funcs]
        l = length(basis.funcs)
        
        # One body integrals
        h = la.Diagonal([E_n(ψ_i) for ψ_i in basis.funcs])
        
        # Two body integrals
        inner = inner_ints(spfs, grid, basis.a)
        u = outer_int(spfs, grid, inner)
        
        # Adding spin
        l = l * 2
        h = kron(h, [1 0; 0 1])
        u = add_spin_u(u)
        spfs = [spfs[(i + 1)÷2] for i in 1:l]
        
        # Anti-symmetrizing u
        u = u .- permutedims(u, [1, 2, 4, 3])
        
        return new{typeof(basis)}(n, l, h, u, spfs, basis, grid)
    end
end

function trapz(f_vals, grid)::Float64 # TODO simpson!!
    val = sum(f_vals)
    val = val - 0.5 * (f_vals[1] + f_vals[end])

    return val * (grid[2] - grid[1])
end

function inner_ints(spfs::Vector{Vector{T}}, grid, a) where T<:Real
    l = length(spfs)
    inner_int = zeros(l, l, length(spfs[1]))
    f_vals = zero(grid)
    coulomb = zero(grid)
    for (xi, x1) in enumerate(grid)
        coulomb .= 1 ./ sqrt.( (grid .- x1).^2 .+ a.^2 )
        for κ in 1:l
            for λ in κ:l
                f_vals .= spfs[κ] .* coulomb .* spfs[λ]
                res = trapz(f_vals, grid)
                inner_int[κ, λ, xi] = res
                inner_int[λ, κ, xi] = res
            end
        end
    end
    return inner_int
end


function inner_ints(spfs::Vector{Vector{T}}, grid, a) where T<:Complex
    l = length(spfs)
    inner_int = zeros(l, l, length(spfs[1]))
    f_vals = zero(grid)
    coloumb = zero(grid)
    for (xi, x1) in enumerate(grid)
        coloumb .= (1 ./ sqrt.( (grid .- x1).^2 .+ a.^2))
        for κ in 1:l
            for λ in κ:l
                f_vals .= conj.(spfs[κ]) .* coloumb .* spfs[λ]
                inner_int[κ, λ, xi] = trapz(f_vals, grid)
            end
        end
    end
    return inner_int
end


function outer_int(spfs, grid, inner_ints)
    l = length(spfs)
    outer_int = zeros(l, l, l, l)
    f_vals = zero(grid)
    
    for κ in 1:l
        for λ in 1:l
            inner = inner_ints[κ, λ, :]
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
            u_new[:, :, ν, λ] .= kron(u_old[:, :, ν_old, λ_old], pattern)
        end
    end
    return u_new
end