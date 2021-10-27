struct RBM2 <: WaveFunction
    M::Int64 # Number of visible units
    N::Int64 # Number of hidden units
    W::Matrix{Float64}
    a::Vector{Float64} # Visible bias
    b::Vector{Float64} # Hidden bias
    σ::Float64
end
function RBM2(visible, hidden, σ, rng::Random.AbstractRNG)
    W = Random.rand(rng, Float64, (visible, hidden))
    a = Random.rand(rng, Float64, visible)
    b = Random.rand(rng, Float64, hidden)
    return RBM2(visible, hidden, W, a, b, σ)
end


function evaluate(particles, wf::RBM2)::Float64
    X = particles.positions
    sum1 = 0.0
    for i in 1:wf.M
        sum1 += (X[i] - wf.a[i])^2
    end
    prod = exp(-sum1 / (2.0 * wf.σ^2))
    for j in 1:wf.N
        sum2 = 0.0
        for i in 1:wf.M
            sum2 += X[i] * wf.W[i, j]
        end
        prod *= 1.0 + exp(wf.b[j] + sum2 / wf.σ^2)
    end
    return sqrt(prod)
end

function kinetic(particles, wf::RBM2)::Float64
    X = particles.positions
    dder = 0.0
    for m in 1:wf.M
        sum1 = 0.0
        sum2 = 0.0
        for n in 1:wf.N
            sum3 = 0.0
            for mm in 1:wf.M
                sum3 += X[mm] * wf.W[mm, n]
            end
            sig_inp = wf.b[n] + sum3 / wf.σ^2
            
            wsig = wf.W[m, n] * sigmoid(sig_inp)
            sum1 += wsig
            sum2 += wsig * wf.W[m, n] * sigmoid(-sig_inp)
        end
        dlnΨ = (wf.a[m] - X[m] + sum1) / 2
        d2lnΨ = (-wf.σ^2 + sum2) / 2
        dder += dlnΨ^2 + d2lnΨ
    end
    return -0.5 * dder / wf.σ^4
end

function QF(particles, idx, wf::RBM2)
    X = particles.positions
    qf = zeros(length(idx)) # maybe use the particles temp vec?
    
    for (i, m) in enumerate(idx)
        sum1 = 0.0
        for n in 1:wf.N
            sum2 = 0.0
            for mm in 1:wf.M
                sum2 += X[mm] * wf.W[mm, n]
            end
            sig_inp = wf.b[n] + sum2 / wf.σ^2
            sum1 += wf.W[m, n] * sigmoid(sig_inp)
        end
        qf[i] = wf.a[m] - X[m] + sum1
    end
    return qf ./ 2
end
function QF!(qf, particles, idx, wf::RBM2)
    X = particles.positions
    for (i, m) in enumerate(idx)
        sum1 = 0.0
        for n in 1:wf.N
            sum2 = 0.0
            for mm in 1:wf.M
                sum2 += X[mm] * wf.W[mm, n]
            end
            sig_inp = wf.b[n] + sum2 / wf.σ^2
            sum1 += wf.W[m, n] * sigmoid(sig_inp)
        end
        qf[i] = wf.a[m] - X[m] + sum1
    end
    return qf ./ 2
end

function paramDer(particles, wf::RBM2)
    X = particles.positions
    da = 1 / wf.σ^2 .* (X .- wf.a) ./ 2
    sig = wf.b
    for m in 1:wf.M
        for n in 1:wf.N
            sig[n] += X[m] .* wf.W[m, n] / wf.σ^2
        end
    end
    sig .= sigmoid.(sig)
    db = sig ./ 2
    dW = X * la.transpose(sig) ./ wf.σ^2 ./ 2
    return da, db, dW
end

function applyGradient(wf::RBM2, grad)
    a = wf.a .- grad[1]
    b = wf.b .- grad[2]
    W = wf.W .- grad[3]
    return RBM(wf.M, wf.N, W, a, b, wf.σ)
end

function sigmoid(x)
    return exp(x) / (1.0 + exp(x))
end