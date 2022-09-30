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

struct ADAM{T} <: Optimizer
    lr::Float64
    β_1::Float64
    m_t::T
    β_2::Float64
    v_t::T
    β_t::Vector{Float64}
    max_iter::Int64
    tol::Float64
    eps::Float64
end

function ADAM(wf; lr, max_iter, tol, β_1=0.9, β_2=0.999)
    m_t = paramDerHolder(wf)
    v_t = paramDerHolder(wf)
    β_t = [β_1, β_2]
    return ADAM(lr, β_1, m_t, β_2, v_t, β_t, max_iter, tol, 1e-8)
end

function update_gradient(gradient, optimizer::ADAM)
    (; lr, β_1, m_t, β_2, v_t, β_t, eps) = optimizer
    β_1t, β_2t = β_t

    for (layer_grad, m_i, v_i) in zip(gradient, m_t, v_t)
        m_i.W_g .*= β_1
        m_i.W_g .+= (1 .- β_1) .* layer_grad.W_g
        m_i.b_g .*= β_1
        m_i.b_g .+= (1 .- β_1) .* layer_grad.b_g

        v_i.W_g .*= β_2
        v_i.W_g .+= (1 .- β_2) .* layer_grad.W_g.^2
        v_i.b_g .*= β_2
        v_i.b_g .+= (1 .- β_2) .* layer_grad.b_g.^2

        layer_grad.W_g .= m_i.W_g ./ (1 .- β_1t) ./ (.√(v_i.W_g ./ (1 .- β_2t)) .+ eps ) .* lr
        layer_grad.b_g .= m_i.b_g ./ (1 .- β_1t) ./ (.√(v_i.b_g ./ (1 .- β_2t)) .+ eps ) .* lr
    end
    β_t[1] *= β_1
    β_t[2] *= β_2

    return gradient
end

function optimize(wf, ham, metro, optimizer; nthreads=1, verbose = true)
    (;max_iter, tol) = optimizer
    E_opt = typemax(Float64)
    wf_opt = wf

    grad_results = Vector{GradientResult}(undef, max_iter) # The sampled values and wavefunctions from each step
    grad_norm = 0
    for i in 1:max_iter
        grad_result = compute_gradient(wf, ham, metro, nthreads) # Running the vmc calculation the get the gradient for this wavefunction
        if grad_result.E < E_opt && grad_result.E > 0
            E_opt = grad_result.E
            wf_opt = wf
        end

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
            print("\rE = $(round(grad_result.E, digits = 6)) iter = $i/$(max_iter)                                      ")
        end
    end
    #if verbose
    #    println("\nNo convergence reached, final norm of gradient was $(grad_norm)")
    #end
    return wf_opt, grad_results
end

function compute_gradient(wf, ham, metro, nthreads)
    distr_steps = distribute_steps(metro.sample_steps, nthreads)
    samplers = [GradientSampler(wf, distr_steps[i]) for i in 1:nthreads]
    
    samplers = steps!(samplers, wf, ham, metro)
    
    grad_result = createResult(samplers)
    return grad_result
end