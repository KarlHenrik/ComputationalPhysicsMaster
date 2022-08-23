abstract type Optimizer end

struct GradientDescent <: Optimizer
    lr::Float64
    max_iter::Int64
    tol::Float64
end

function update_gradient(gradient, optimizer::GradientDescent)
    return gradient .* optimizer.lr
end

function optimize(wf, ham, metro, optimizer; threads)
    (;max_iter, tol) = optimizer
    
    grad_results = Vector{GradientResult}(undef, maxiter) # The sampled values and wavefunctions from each step
    grad_norm = 0
    for i in 1:maxiter
        grad_result = compute_gradient(wf, ham, metro) # Running the vmc calculation the get the gradient for this wavefunction
        grad_results[i] = grad_result
        
        grad = update_gradient!(optimizer, grad_result.gradient)
        
        grad_norm = la.norm(grad)
        if grad_norm < tol
            grad_results = grad_results[1:i]
            return wf, grad_results
        end
        
        wf = applyGradient(wf, grad)
        print("\rE = $(round(grad_result.E, digits = 3)) iter = $i/$(maxiter)                       ")
    end
    println("No convergence reached, final norm of gradient was $(grad_norm)")
    return wf, grad_results
end

function compute_gradient(wf, ham, metro, nthreads)
    distr_steps = distribute_steps(metro.sampled_steps, nthreads)
    samplers = [GradientSampler(wf, distr_steps[i]) for i in 1:nthreads]
    
    samplers = vmc!(samplers, wf, ham, metro)
    
    grad_result = CreateResult(samplers)
    return grad_result
end

struct GradientResult{T} <: Result
    E::Float64
    E_std::Float64
    E2::Float64
    gradient::T
end

function createResult(samplers::Vector{GradientSampler{T}}) where T
    E = 0.0
    E2 = 0.0
    ∇Ψ = zero.(samplers[1].∇Ψ)
    ∇ΨE = zero.(samplers[1].∇ΨE)
    
    totalSteps = 0
    for sampler in samplers
        totalSteps += sampler.sampled_steps
    end
    
    for sampler in samplers
        E += sampler.E / totalSteps
        E2 += sampler.E2 / totalSteps
        ∇Ψ = ∇Ψ .+ sampler.∇Ψ ./ totalSteps
        ∇ΨE = ∇ΨE .+ sampler.∇ΨE ./ totalSteps
    end
    E_std = sqrt(abs(E2 - E^2))
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end

function createResult(sampler::GradientSampler)
    E = sampler.E / sampler.sampled_steps
    E2 = sampler.E2 / sampler.sampled_steps
    ∇Ψ = sampler.∇Ψ ./ sampler.sampled_steps
    ∇ΨE = sampler.∇ΨE ./ sampler.sampled_steps
    
    E_std = sqrt(E2 - E^2)
    gradient = 2 .* (∇ΨE .- E .* ∇Ψ)
    
    return GradientResult(E, E_std, E2, gradient)
end