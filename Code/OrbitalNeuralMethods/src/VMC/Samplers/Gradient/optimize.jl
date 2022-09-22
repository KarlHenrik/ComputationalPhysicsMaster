abstract type Optimizer end

struct GradientDescent <: Optimizer
    lr::Float64
    max_iter::Int64
    tol::Float64
    function GradientDescent(;lr, max_iter, tol)
        return new(lr, max_iter, tol)
    end
end

function update_gradient(gradient, optimizer::GradientDescent)
    return gradient .* optimizer.lr
end

function optimize(wf, ham, metro, optimizer; nthreads=1, verbose = true)
    (;max_iter, tol) = optimizer
    
    grad_results = Vector{GradientResult}(undef, max_iter) # The sampled values and wavefunctions from each step
    grad_norm = 0
    for i in 1:max_iter
        grad_result = compute_gradient(wf, ham, metro, nthreads) # Running the vmc calculation the get the gradient for this wavefunction
        grad_results[i] = grad_result
        
        grad = update_gradient(grad_result.gradient, optimizer)
        
        grad_norm = la.norm(grad)
        if grad_norm < tol || isnan(grad_result.E)
            grad_results = grad_results[1:i]
            println("\nConvergence reached, final norm of gradient was $(grad_norm)")
            return wf, grad_results
        end
        
        wf = applyGradient(wf, grad)
        if verbose
            print("\rE = $(round(grad_result.E, digits = 3)) iter = $i/$(max_iter)                                      ")
        end
    end
    if verbose
        println("\nNo convergence reached, final norm of gradient was $(grad_norm)")
    end
    
    return wf, grad_results
end

function compute_gradient(wf, ham, metro, nthreads)
    distr_steps = distribute_steps(metro.sample_steps, nthreads)
    samplers = [GradientSampler(wf, distr_steps[i]) for i in 1:nthreads]
    
    samplers = steps!(samplers, wf, ham, metro)
    
    grad_result = createResult(samplers)
    return grad_result
end