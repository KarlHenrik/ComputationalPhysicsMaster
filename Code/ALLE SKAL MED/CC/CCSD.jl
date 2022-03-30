struct CCSDState{T}
    system::T
    ϵ::Vector{Float64}
    
    α::Float64
    
    t1::Matrix{Float64}
    t1_new::Matrix{Float64}
    
    t2::Array{Float64, 4}
    t2_new::Array{Float64, 4}
end

function setup_CCSD(system; α)
    (; l, h, u) = system
    
    ϵ = zeros(l)
    @inbounds for q in 1:l
        ϵ[q] = h[q, q]
        for i in 1:n
            ϵ[q] += u[q, i, q, i]
        end
    end
    
    t1_new = zeros((l, l))
    t1 = zeros((l, l))
    @inbounds for i in 1:n
        for a in n+1:l
            t1[a, i] = h[a, i] / (ϵ[i] - ϵ[a])
        end
    end
    
    t2_new = zeros((l, l, l, l))
    t2 = zeros((l, l, l, l))
    @inbounds for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    t2[a, b, i, j] = u[a, b, i, j] / (ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b])
                end
            end
        end
    end
    
    return CCSDState(system, ϵ, α, t1, t1_new, t2, t2_new)
end

function energy(state::CCSDState)
    #=
    Returns E_CCSD - E_0
    
    where E_0 is the energy of the reference determinant
    =#
    
    (; system, t1, t2) = state
    (; n, l, h, u) = system
    
    E = 0.0
    
    @inbounds for i in 1:n
        for a in n+1:l
            E += h[i, a] * t1[a, i]
        end
    end
    
    @inbounds for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    E += 0.25 * u[i, j, a, b] * t2[a, b, i, j]
                    E += 0.5 * u[i, j, a, b] * t1[a, i] * t1[b, j]
                end
            end
        end
    end
    
    return E
end



function CCSD_Update!(state::CCSDState)
    (; system, α, ϵ, t1, t1_new, t2, t2_new) = state
    (; n, l, h, u) = system
    
    L = system.l # There is a collision of notation. In here, L is the number of basis functions
    
    """
    The t1 amplitude equations, from page 75 of Crawford & Schaefer
    """
    
    @inbounds Threads.@threads for a in n+1:L
        for i in 1:n
            x = 0.0

            x += h[a, i]

            # Terms that were moved over
            x -= h[a, a] * t1[a, i]
            x += h[i, i] * t1[a, i]

            for c in n+1:L
                x += h[a, c] * t1[c, i]
            end

            for k in 1:n
                x -= h[k, i] * t1[a, k]
            end

            for k in 1:n
                for c in n+1:L
                    x += u[k, a, c, i] * t1[c, k]

                    x += h[k, c] * t2[a, c, i, k]

                    x -= h[k, c] * t1[c, i] * t1[a, k]
                end
            end

            for k in 1:n
                for c in n+1:L
                    for d in n+1:L
                        x += 0.5 * u[k, a, c, d] * t2[c, d, k, i]

                        x -= u[k, a, c, d] * t1[c, k] * t1[d, i]
                    end
                end
            end

            for k in 1:n
                for l in 1:n
                    for c in n+1:L
                        x -= 0.5 * u[k, l, c, i] * t2[c, a, k, l]

                        x -= u[k, l, c, i] * t1[c, k] * t1[a, l]
                    end
                end
            end

            for k in 1:n
                for l in 1:n
                    for c in n+1:L
                        for d in n+1:L
                            x -= u[k, l, c, d] * t1[c, k] * t1[d, i] * t1[a, l]

                            x += u[k, l, c, d] * t1[c, k] * t2[d, a, l, i]

                            x -= 0.5 * u[k, l, c, d] * t2[c, d, k, i] * t1[a, l]

                            x -= 0.5 * u[k, l, c, d] * t2[c, a, k, l] * t1[d, i]
                        end
                    end
                end
            end

            ϵ_ai = ϵ[i] - ϵ[a]
            x = α * t1[a, i] + (1 - α) * x / ϵ_ai
            t1_new[a, i] = x
        end
    end
    
    """
    The t2 amplitude equations, from page 76 of Crawford & Schaefer
    """
    @inbounds Threads.@threads for a in n+1:L
        for i in 1:n
            for j in i+1:n
                for b in a+1:L
                    s = 0.0

                    s += u[a, b, i, j]
                    
                    # Terms that were moved over
                    s -= h[b, b] * t2[a, b, i, j] # -δbc
                    s += h[a, a] * t2[b, a, i, j] # Pab δbc
                    s += h[j, j] * t2[a, b, i, j] # δkj
                    s -= h[i, i] * t2[a, b, j, i] # -Pij δkj

                    for c in n+1:L
                        s += h[b, c] * t2[a, c, i, j] # 1
                        s -= h[a, c] * t2[b, c, i, j] # -Pab
                        
                        s += u[a, b, c, j] * t1[c, i] # 1
                        s -= u[a, b, c, i] * t1[c, j] # -Pij
                    end
                    
                    for k in 1:n
                        s -= h[k, j] * t2[a, b, i, k] # -1
                        s += h[k, i] * t2[a, b, j, k] # Pij
                        
                        s -= u[k, b, i, j] * t1[a, k] # -1
                        s += u[k, a, i, j] * t1[b, k] # Pab
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            s += 0.5 * u[k, l, i, j] * t2[a, b, k, l]
                            
                            s += 0.5 * u[k, l, i, j] * t1[a, k] * t1[b, l] # 1
                            s -= 0.5 * u[k, l, i, j] * t1[b, k] * t1[a, l] # -Pab
                        end
                    end
                    
                    for c in n+1:L
                        for d in n+1:L
                            s += 0.5 * u[a, b, c, d] * t2[c, d, i, j]
                            
                            s += 0.5 * u[a, b, c, d] * t1[c, i] * t1[d, j] # 1
                            s -= 0.5 * u[a, b, c, d] * t1[c, j] * t1[d, i] # -Pij
                        end
                    end
                    
                    for k in 1:n
                        for c in n+1:L
                            s += u[k, b, c, j] * t2[a, c, i, k] # 1
                            s -= u[k, b, c, i] * t2[a, c, j, k] # -Pij
                            s -= u[k, a, c, j] * t2[b, c, i, k] # -Pab
                            s += u[k, a, c, i] * t2[b, c, j, k] # Pij Pab
                            
                            s -= u[k, b, i, c] * t1[a, k] * t1[c, j] # -1
                            s += u[k, b, j, c] * t1[a, k] * t1[c, i] # Pij
                            s += u[k, a, i, c] * t1[b, k] * t1[c, j] # Pab
                            s -= u[k, a, j, c] * t1[b, k] * t1[c, i] # -Pij Pab
                            
                            s += h[k, c] * t1[a, k] * t2[b, c, i, j] # 1
                            s -= h[k, c] * t1[b, k] * t2[a, c, i, j] # -Pab
                            
                            s += h[k, c] * t1[c, i] * t2[a, b, j, k] # 1
                            s -= h[k, c] * t1[c, j] * t2[a, b, i, k] # -Pij
                        end
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            for c in n+1:L
                                s -= u[k, l, c, i] * t1[c, k] * t2[a, b, l, j] # -1
                                s += u[k, l, c, j] * t1[c, k] * t2[a, b, l, i] # Pij
                                
                                s += u[k, l, i, c] * t1[a, l] * t2[b, c, j, k] # 1
                                s -= u[k, l, j, c] * t1[a, l] * t2[b, c, i, k] # -Pij
                                s -= u[k, l, i, c] * t1[b, l] * t2[a, c, j, k] # -Pab
                                s += u[k, l, j, c] * t1[b, l] * t2[a, c, i, k] # Pij Pab
                                
                                s += 0.5 * u[k, l, c, j] * t1[c, i] * t2[a, b, k, l] # 1
                                s -= 0.5 * u[k, l, c, i] * t1[c, j] * t2[a, b, k, l] # -Pij
                                
                                s += 0.5 * u[k, l, c, j] * t1[c, i] * t1[a, k] * t1[b, l] # 1
                                s -= 0.5 * u[k, l, c, i] * t1[c, j] * t1[a, k] * t1[b, l] # -Pij
                                s -= 0.5 * u[k, l, c, j] * t1[c, i] * t1[b, k] * t1[a, l] # -Pab
                                s += 0.5 * u[k, l, c, i] * t1[c, j] * t1[b, k] * t1[a, l] # Pij Pab
                            end
                        end
                    end
                    
                    for k in 1:n
                        for c in n+1:L
                            for d in n+1:L
                                s += u[k, a, c, d] * t1[c, k] * t2[d, b, i, j] # 1
                                s -= u[k, b, c, d] * t1[c, k] * t2[d, a, i, j] # -Pab
                                
                                s += u[a, k, d, c] * t1[d, i] * t2[b, c, j, k] # 1
                                s -= u[a, k, d, c] * t1[d, j] * t2[b, c, i, k] # -Pij
                                s -= u[b, k, d, c] * t1[d, i] * t2[a, c, j, k] # -Pab
                                s += u[b, k, d, c] * t1[d, j] * t2[a, c, i, k] # Pij Pab
                                
                                s -= 0.5 * u[k, b, c, d] * t1[a, k] * t2[c, d, i, j] # -1
                                s += 0.5 * u[k, a, c, d] * t1[b, k] * t2[c, d, i, j] # Pab
                                
                                s -= u[k, b, c, d] * t1[c, i] * t1[a, k] * t1[d, j] # -1
                                s += u[k, b, c, d] * t1[c, j] * t1[a, k] * t1[d, i] # Pij
                                s += u[k, a, c, d] * t1[c, i] * t1[b, k] * t1[d, j] # Pab
                                s -= u[k, a, c, d] * t1[c, j] * t1[b, k] * t1[d, i] # -Pij Pab
                            end
                        end
                    end
                    
                    for k in 1:n
                        for l in 1:n
                            for c in n+1:L
                                for d in n+1:L
                                    s += 0.5 * u[k, l, c, d] * t2[a, c, i, k] * t2[d, b, l, j] # 1
                                    s -= 0.5 * u[k, l, c, d] * t2[a, c, j, k] * t2[d, b, l, i] # -Pij
                                    s -= 0.5 * u[k, l, c, d] * t2[b, c, i, k] * t2[d, a, l, j] # -Pab
                                    s += 0.5 * u[k, l, c, d] * t2[b, c, j, k] * t2[d, a, l, i] # Pij Pab
                                    
                                    s += 0.25 * u[k, l, c, d] * t2[c, d, i, j] * t2[a, b, k, l]

                                    s -= 0.5 * u[k, l, c, d] * t2[a, c, i, j] * t2[b, d, k, l] # -1
                                    s += 0.5 * u[k, l, c, d] * t2[b, c, i, j] * t2[a, d, k, l] # Pab
                                    
                                    s -= 0.5 * u[k, l, c, d] * t2[a, b, i, k] * t2[c, d, j, l] # -1
                                    s += 0.5 * u[k, l, c, d] * t2[a, b, j, k] * t2[c, d, i, l] # Pij
                                    
                                    s -= u[k, l, c, d] * t1[c, k] * t1[d, i] * t2[a, b, l, j] # -1
                                    s += u[k, l, c, d] * t1[c, k] * t1[d, j] * t2[a, b, l, i] # Pij
                                    
                                    s -= u[k, l, c, d] * t1[c, k] * t1[a, l] * t2[d, b, i, j] # -1
                                    s += u[k, l, c, d] * t1[c, k] * t1[b, l] * t2[d, a, i, j] # Pab
                                    
                                    s += 0.25 * u[k, l, c, d] * t1[c, i] * t1[d, j] * t2[a, b, k, l] # 1
                                    s -= 0.25 * u[k, l, c, d] * t1[c, j] * t1[d, i] * t2[a, b, k, l] # -Pij
                                    
                                    s += 0.25 * u[k, l, c, d] * t1[a, k] * t1[b, l] * t2[c, d, i, j] # 1
                                    s -= 0.25 * u[k, l, c, d] * t1[b, k] * t1[a, l] * t2[c, d, i, j] # -Pab
                                    
                                    s += u[k, l, c, d] * t1[c, i] * t1[b, l] * t2[a, d, k, j] # 1
                                    s -= u[k, l, c, d] * t1[c, j] * t1[b, l] * t2[a, d, k, i] # -Pij
                                    s -= u[k, l, c, d] * t1[c, i] * t1[a, l] * t2[b, d, k, j] # -Pab
                                    s += u[k, l, c, d] * t1[c, j] * t1[a, l] * t2[b, d, k, i] # Pij Pab
                                    
                                    s += 0.25 * u[k, l, c, d] * t1[c, i] * t1[a, k] * t1[d, j] * t1[b, l] # 1
                                    s -= 0.25 * u[k, l, c, d] * t1[c, j] * t1[a, k] * t1[d, i] * t1[b, l] # -Pij
                                    s -= 0.25 * u[k, l, c, d] * t1[c, i] * t1[b, k] * t1[d, j] * t1[a, l] # -Pab
                                    s += 0.25 * u[k, l, c, d] * t1[c, j] * t1[b, k] * t1[d, i] * t1[a, l] # Pij Pab
                                end
                            end
                        end
                    end
                    
                    ϵ_abij = ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b]
                    s = α * t2[a, b, i, j] + (1 - α) * s / ϵ_abij
                    t2_new[a, b, i, j] = s
                    t2_new[a, b, j, i] = -s
                    t2_new[b, a, i, j] = -s
                    t2_new[b, a, j, i] = s
                end
            end
        end
    end
    
    t1 .= t1_new
    t2 .= t2_new
    
    return state
end
;