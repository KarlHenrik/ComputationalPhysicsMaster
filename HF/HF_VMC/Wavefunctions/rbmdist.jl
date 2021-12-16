struct RBMDist{T} <: WaveFunction
    # Normal rbm stuff
    M::Int64 # Number of visible units
    N::Int64 # Number of hidden units
    W::Matrix{Float64}
    a::Vector{Float64} # Visible bias
    b::Vector{Float64} # Hidden bias
    σ::Float64
    
    # Hessian setup
    hes_grad_result::Vector{Float64}
    hes_grad_tapes::Dict{DataType, Any}
    hes_jac_result::Matrix{Float64}
    hes_jac_config::T
end
function RBMDist(visible, hidden, σ, rng::Random.AbstractRNG)
    W = Random.rand(rng, Float64, (visible, hidden))
    a = Random.rand(rng, Float64, visible)
    b = Random.rand(rng, Float64, hidden)

    n = visible - 1

    hes_grad_result = zeros(n)
    config = rd.GradientConfig(hes_grad_result)
    hes_grad_tapes = Dict{DataType, Any}()
    hes_jac_result = zeros(n, n)
    hes_jac_config = fwdd.JacobianConfig(nothing, hes_grad_result, hes_grad_result)

    return RBMDist{typeof(hes_jac_config)}(visible, hidden, W, a, b, σ,
                                       hes_grad_result, hes_grad_tapes, hes_jac_result, hes_jac_config)
end

function RBMDist(wf::RBMDist, W, a, b)
    W = W
    a = a
    b = b
    return RBMDist{typeof(wf.hes_jac_config)}(wf.M, wf.N, W, a, b, wf.σ,
                                           wf.hes_grad_result, wf.hes_grad_tapes, wf.hes_jac_result, wf.hes_jac_config)
end

function evaluate(particles, wf::RBMDist)
    return evaluate(particles.positions, wf)
end

function evaluate(X::AbstractArray, wf::RBMDist)
    r = dist(X)
    sum1 = 0.0
    for i in 1:wf.M-1
        sum1 += (X[i] - wf.a[i])^2
    end
    sum1 += (r - wf.a[end])^2
    prod = exp(-sum1 / (2.0 * wf.σ^2))
    for j in 1:wf.N
        sum2 = 0.0
        for i in 1:wf.M-1
            sum2 += X[i] * wf.W[i, j]
        end
        sum2 += r * wf.W[end, j]
        prod *= 1.0 + exp(wf.b[j] + sum2 / wf.σ^2)
    end
    return prod
end

function hes_grad!(y, x::Array{T}, wf::RBMDist) where {T<:Real}
    if !haskey(wf.hes_grad_tapes, T)
        config = rd.GradientConfig(x)
        tape = rd.compile(rd.GradientTape(x -> evaluate(x, wf), x, config))
        wf.hes_grad_tapes[T] = tape
    end
    tape = wf.hes_grad_tapes[T]
    return rd.gradient!(y, tape, x)
end

function kinetic(particles, wf::RBMDist)
    x = particles.positions
    fwdd.jacobian!(wf.hes_jac_result, (y, x) -> hes_grad!(y, x, wf), wf.hes_grad_result, x, wf.hes_jac_config)
    return -0.5 * la.tr(wf.hes_jac_result) / evaluate(x, wf)
end

function paramDer(particles, wf::RBMDist)
    X = vcat(particles.positions, dist(particles.positions))
    da = 1 / wf.σ^2 * (X .- wf.a)
    sig = wf.b
    for m in 1:wf.M
        for n in 1:wf.N
            sig[n] += X[m] .* wf.W[m, n] / wf.σ^2
        end
    end
    sig .= sigmoid.(sig)
    db = sig
    dW = X * la.transpose(sig) ./ wf.σ^2
    return da, db, dW
end

function applyGradient(wf::RBMDist, grad)
    a = wf.a .- grad[1]
    b = wf.b .- grad[2]
    W = wf.W .- grad[3]
    return RBMDist(wf, W, a, b)
end

function sigmoid(x)
    return exp(x) / (1.0 + exp(x))
end

function dist(p)
    return (p[1] - p[3])^2 + (p[2] - p[4])^2
end