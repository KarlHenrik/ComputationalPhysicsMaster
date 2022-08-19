function compute_ground_state!(state; max_iter = 1000, tol = 1e-10, verbose = false)
    old_E = energy(state)
    delta_E = 0
    for i in 1:max_iter
        update!(state)
        
        new_E = energy(state)
        delta_E = old_E - new_E
        if verbose
            println("Iteration $i: E = $new_E")
        end
        if abs(delta_E) < tol
            return state
        end
        old_E = new_E
    end
    println("Did not converge after $(max_iter) iterations. Final energy change was $(delta_E)")
    return state
end