struct HFState{T}
    system::T
    
    C::Matrix{Float64}
    P::Matrix{Float64}
    F::Matrix{Float64}
end

function HF(system)
    (; l) = system
    
    C = la.I(l)
    P = zeros((l, l))
    F = zeros((l, l))
    
    state = HFState{typeof(system)}(system, C, P, F)
    P_update!(state)
    F_update!(state)
    return state
end

function HF_update!(state::HFState; iters = 1)
    (; C, F) = state
    for i in 1:iters
        P_update!(state)
        F_update!(state)
        C .= la.eigvecs(F)
    end
    return state
end

function update!(state::HFState)
    (; C, F) = state
    
    C .= la.eigvecs(F)
    P_update!(state)
    F_update!(state)
        
    return state
end

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
    (; l, h, u) = state.system
    
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
;