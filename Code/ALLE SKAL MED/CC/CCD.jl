struct CCD_System
    n::Int64
    l::Int64
    α::Float64
    
    h::Matrix{Float64}
    u::Array{Float64, 4}
    ϵ::Vector{Float64}
    
    t::Array{Float64, 4}
    t_new::Array{Float64, 4}
end

function setup_CCD(h, u, n; α)
    l = size(h)[1]
    
    ϵ = zeros(l)
    for q in 1:l
        ϵ[q] = h[q, q]
        for i in 1:n
            ϵ[q] += u[q, i, q, i]
        end
    end
    
    t = zeros((l, l, l, l))
    t_new = zeros((l, l, l, l))
    
    for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    t[a, b, i, j] = u[a, b, i, j] / (ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b])
                end
            end
        end
    end
    
    return CCD_System(n, l, α, h, u, ϵ, t, t_new)
end

function E_CCD(system)
    (; n, l, u, t) = system
    
    E = 0.0
    for i in 1:n
        for j in 1:n
            for a in n+1:l
                for b in n+1:l
                    E += u[i, j, a, b] * t[a, b, i, j]
                end
            end
        end
    end
    return 0.25 * E
end

function CCD_Update!(system)
    (; n, l, u, ϵ, t, t_new, α) = system
    
    for i in 1:n
        for j in i+1:n
            for a in n+1:l
                for b in a+1:l
                    # Page 361 of An Advanced Course in Computational Physics
                    s = 0

                    s += u[a, b, i, j]

                    for c in n+1:l
                        for d in n+1:l
                            s += 0.5 * u[a, b, c, d] * t[c, d, i, j]
                        end
                    end

                    for k in 1:n
                        for li in 1:n
                            s += 0.5 * u[k, li, i, j] * t[a, b, k, li]
                        end
                    end

                    for k in 1:n
                        for c in n+1:l
                            s += u[k, b, c, j] * t[a, c, i, k] # 1
                            s -= u[k, b, c, i] * t[a, c, j, k] # -Pij
                            s -= u[k, a, c, j] * t[b, c, i, k] # -Pab
                            s += u[k, a, c, i] * t[b, c, j, k] # Pij Pab
                        end
                    end

                    for k in 1:n
                        for li in 1:n
                            for c in n+1:l
                                for d in n+1:l
                                    s += 0.25 * u[k, li, c, d] * t[c, d, i, j] * t[a, b, k, li]

                                    s += u[k, li, c, d] * t[a, c, i, k] * t[b, d, j, li] # 1
                                    s -= u[k, li, c, d] * t[a, c, j, k] * t[b, d, i, li] # -Pij

                                    s -= 0.5 * u[k, li, c, d] * t[d, c, i, k] * t[a, b, li, j] # -1
                                    s += 0.5 * u[k, li, c, d] * t[d, c, j, k] * t[a, b, li, i] # Pij

                                    s -= 0.5 * u[k, li, c, d] * t[a, c, li, k] * t[d, b, i, j] # -1
                                    s += 0.5 * u[k, li, c, d] * t[b, c, li, k] * t[d, a, i, j] # Pab
                                end
                            end
                        end
                    end

                    ϵ_abij = ϵ[i] + ϵ[j] - ϵ[a] - ϵ[b]
                    s = α * t[a, b, i, j] + (1 - α) * s / ϵ_abij
                    t_new[a, b, i, j] = s
                    t_new[a, b, j, i] = -s
                    t_new[b, a, i, j] = -s
                    t_new[b, a, j, i] = s
                end
            end
        end
    end
    t .= t_new
    return system
end
;