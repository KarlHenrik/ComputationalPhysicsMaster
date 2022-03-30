struct Slater{T}
    C::Matrix{Float64}
    l::Int64
    basis::T
end
function Slater(state::HFState)
    (; C, system) = state
    (; l, basis) = system
    
    return Slater(C, l, basis)
end

function evaluate!(walker, x, wf::Slater{Basis})
    (; basis_eval) = walker.wf_mut
    basis_eval .= evaluate!(basis_eval, x, wf.basis) # the basis functions evaluated at x
end

function evaluate!(walker, x, wf::Slater{SpinBasis})
    (; basis_eval, nospin) = walker.wf_mut
    
    nospin .= evaluate!(nospin, x, wf.basis.base) # the basis functions evaluated at x
    @inbounds for i in eachindex(nospin) # doubling the basis functions to include spin
        basis_eval[2i-1] = nospin[i]
        basis_eval[2i] = nospin[i]
    end
end

function consider!(walker, system, p_idx, dim, dist)
    wf = system.wf
    x = walker.positions[p_idx] + dist
    basis_eval = evaluate!(walker, x, wf) # the basis functions evaluated at x
    
    (; slater, old_slater_col) = walker.wf_mut
    (; C, l) = wf
    for col in 1:system.n
        ϕ_col_x = 0.0
        for j in 1:l
            ϕ_col_x += C[j, col] * basis_eval[j]
        end
        old_slater_col[i] = slater[i, col]
        slater[i, col] = ϕ_col_x
    end
    
    return la.det(slater)
end