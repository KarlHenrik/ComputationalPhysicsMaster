using LinearAlgebra: norm

# Schemes for optimizing the wavefunction
abstract type OptimizerScheme <: Scheme end

createSampler(scheme::OptimizerScheme, wf, sampled_steps) = GradientSampler(wf, sampled_steps)

struct GradientDescent <: OptimizerScheme
    lr::Float64
    maxiter::Int64
    tol::Float64
end
GradientDescent(;lr, maxiter, tol) = GradientDescent(lr, maxiter, tol)

function update_gradient(gradient, scheme::GradientDescent)
    return gradient .* scheme.lr
end

function run_scheme(wf, ham, metro, dims, num, nthreads, scheme::OptimizerScheme)
    # Saving the sampled values and wavefunctions from each step to analyze later
    results = Vector{GradientResult}(undef, scheme.maxiter)
    wavefunctions = Vector{WaveFunction}(undef, scheme.maxiter)

    for i in 1:scheme.maxiter
        # Running the vmc calculation the get the gradient for this wavefunction
        result = vmc(wf, ham, metro, dims, num, nthreads, scheme)
        # Saving the results for this step for later analysis
        results[i] = result
        wavefunctions[i] = wf
        # Computing the gradient using a chosen optimization method
        grad = update_gradient(result.gradient, scheme) # This is where the optimization scheme is applied (Gradient Descent, Adam, etc.)
        # Stops the calculation if the gradient is smaller than the given tolerance
        if norm(grad) < scheme.tol
            results = results[1:i]
            wavefunctions = wavefunctions[1:i]
            break
        end
        # Applying gradient
        wf = applyGradient(wf, grad)

        # Output to keep track of progress
        if i % 5 == 0
            print("\rE = $(round(samples.E, digits = 3)) iter = $i/$(opt.maxiter)")
        end
    end
    
    return results, wavefunctions
end