function getP(C, system)
    (; n, l) = system
    
    P = la.zeros((l, l))
    for i in 1:n
        for a in 1:l
            for b in 1:l
                @inbounds P[b, a] += conj(C[a, i]) * C[b, i]
            end
        end
    end
    return P
end

function getF(P, system)
    (; h, u, n, l) = system
    
    F = la.zeros((l, l))
    F .+= h
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

function SCF(C, system, iters)
    for i in 1:iters
        P = getP(C, system)
        F = getF(P, system)
        C = la.eigvecs(F)
    end
    
    return C
end