struct RHFState{T}
    system::SpatialSystem{SpinBasis{T}}
    
    # Restricted system
    l::Int64
    h::Matrix{Float64} # compact one-body integrals
    u::Array{Float64, 4} # compact two-body integrals
    
    # Hartree-Fock related coefficients
    C::Matrix{Float64}
    P::Matrix{Float64}
    F::Matrix{Float64}
end

function setup_RHF(system::SpatialSystem{SpinBasis{T}}) where T <: SpatialBasis
    (; transform, basis, h, grid, V) = system
    
    @assert system.n % 2 == 0 "Closed shell restricted systems only accept an even number of electrons"
    @assert la.I(size(transform)[1]) == transform "Cannot use transformed system in RHF"
    
    l = basis.base.l
    h = shrink_restricted(h)
    spfs = spatial(basis.base, grid)
    inner = inner_ints(spfs, grid, V)
    u = outer_int(spfs, grid, inner)
    u = 2 .* u .- permutedims(u, [1, 2, 4, 3])
    
    C = la.I(l)
    P = zeros((l, l))
    F = zeros((l, l))
    
    state = RHFState{typeof(system.basis.base)}(system, l, h, u, C, P, F)
    P_update!(state)
    F_update!(state)
    return state
end

function update!(state::RHFState)
    (; C, F) = state
    
    C .= la.eigvecs(F)
    P_update!(state)
    F_update!(state)
        
    return state
end

# Szabo p.139
function P_update!(state::RHFState)
    (; P, C, l) = state
    (; n) = state.system
    
    for a in 1:l
        for b in 1:l
            @inbounds P[b, a] = 0
        end
    end
    
    for i in 1:n÷2
        for a in 1:l
            for b in 1:l
                @inbounds P[b, a] += 2 * conj(C[a, i]) * C[b, i]
            end
        end
    end
    return P
end

# Szabo p.141
function F_update!(state::RHFState)
    (; P, F, l, h, u) = state
    
    F .= h
    for c in 1:l
        for d in 1:l
            @inbounds P_dc = P[d, c]
            for a in 1:l
                for b in 1:l
                    @inbounds F[a, b] += P_dc * 0.5 * u[a, c, b, d]
                end
            end
        end
    end
    return F
end

# Szabo p.150
function energy(state::RHFState)
    (; P, F, l ,h, u) = state
    
    energy = 0.0
    for a in 1:l
        for b in 1:l
            @inbounds energy += 0.5 * P[b, a] * (h[a, b] + F[a, b])
        end
    end
    return real(energy)
end

function System(state::RHFState)
    return System(state.system, expand_restricted(state.C))
end

#=
[A_11 0    A_12 0
 0    A_11 0    A_12
 A_21 0    A_22 0
 0    A_21 0    A_22]
=>
[A_11 A_12
 A_21 A_22]
=#
function shrink_restricted(A)
    l = size(A)[1] ÷ 2
    
    S = zeros(l, l)
    for i in 1:l
        for j in 1:l
            S[i, j] = A[i*2-1, j*2-1]
        end
    end
    return S
end

#=
[A_11 A_12
 A_21 A_22]
=>
[A_11 0    A_12 0
 0    A_11 0    A_12
 A_21 0    A_22 0
 0    A_21 0    A_22]
=#
function expand_restricted(A)
    l = size(A)[1]
    
    S = zeros(2 * l, 2 * l)
    for i in 1:l
        for j in 1:l
            S[2*i-1, 2*j-1] = A[i, j]
            S[2*i, 2*j] = A[i, j]
        end
    end
    return S
end

function checkRestricted(C)
    l = size(C)[1]
    if l%2 != 0
        return false
    end
    # Only nonzero values at even rows at odd columns and odd columns at even rows
    for i in 2:2:l
        for j in 1:2:l
            if C[i, j] != 0
                return false
            end
            if C[j, i] != 0
                return false
            end
        end
    end
    # Every value reappears at its index +1 +1
    for i in 1:2:l
        for j in 1:2:l
            if C[i, j] != C[i+1, j+1]
                return false
            end
        end
    end

    return true
end

#= Not in use
function check_restricted(C)
    l = size(C)[1]
    
    restricted = true
    for func in 1:2:l
        for coef in 1:2:l
            orthogonal = (C[coef, func] == C[coef+1, func+1]) # Up coefficient == down coefficient
            odd_up = C[coef+1, func] == 0 # Function 1 has no down part
            even_down = C[coef, func+1] == 0 # Function 2 has no up part
            if !(orthogonal & odd_up & even_down)
                restricted = false
            end
        end
    end
    return restricted
end
=#