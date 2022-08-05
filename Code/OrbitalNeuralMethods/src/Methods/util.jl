function compute_ground_state!(state; max_iter = 1000, tol = 1e-10, verbose = false)
    old_E = energy(state)
    
    for i in 1:max_iter
        update!(state)
        
        new_E = energy(state)
        if verbose
            println("Iteration $i: E = $new_E")
        end
        if abs(new_E - old_E) < tol
            return state
        end
        old_E = new_E
    end
    return state
end