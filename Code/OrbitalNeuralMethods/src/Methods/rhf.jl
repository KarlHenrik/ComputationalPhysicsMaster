struct RHFState{T}
    system::SpatialSystem{SpinBasis{T}}
    
    h_r::Matrix{Float64} # compact one-body integrals
    u_r::Array{Float64, 4} # compact two-body integrals
    C::Matrix{Float64}
    P::Matrix{Float64}
    F::Matrix{Float64}
end

function setup_RHF(system::SpatialSystem{SpinBasis{T}})
    (; transform, basis, h, grid, V) = system
    
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
    
    state = RHFState{typeof(system)}(system, h, u, C, P, F)
    P_update!(state)
    F_update!(state)
end

function HF_update!(state::HFState; iters)
    (; C, F) = state
    for i in 1:iters
        P_update!(state)
        F_update!(state)
        C .= la.eigvecs(F)
    end
    return state
end
HF_update!(state::HFState) = HF_update!(state, iters = 1)

function P_update!(state::HFState)
    (; P, C) = state
    (; n, l) = state.system
    
    for a in 1:l
        for b in 1:l
            @inbounds P[b, a] = 0
        end
    end
    
    for i in 1:n
        for a in 1:l
            for b in 1:l
                @inbounds P[b, a] += conj(C[a, i]) * C[b, i]
            end
        end
    end
    return P
end

function F_update!(state::HFState)
    (; P, F) = state
    (; n, l, h, u) = state.system
    
    F .= h
    for c in 1:l
        for d in 1:l
            @inbounds P_dc = P[d, c]
            for a in 1:l
                for b in 1:l
                    @inbounds F[a, b] += P_dc * u[a, c, b, d]
                end
            end
        end
    end
    return F
end

function energy(state::HFState)
    (; P) = state
    (; l, h, u) = state.system
    
    energy = 0.0
    for a in 1:l
        for b in 1:l
            @inbounds energy += P[b, a] * h[a, b]
            for c in 1:l
                for d in 1:l
                    @inbounds energy += 0.5 * P[b, a] * P[d, c] * u[a, c, b, d]
                end
            end
        end
    end
    return real(energy)
end

function System(state::HFState)
    return System(state.system, state.C)
end

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