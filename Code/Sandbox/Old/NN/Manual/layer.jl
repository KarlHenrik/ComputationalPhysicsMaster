abstract type Layer end

function layerParams(layer::Layer)
    return -1, -1
end

# --------------------- Dense -----------------------
struct Dense <: Layer
    W::Matrix{Float64}
    b::Vector{Float64}
end
Dense(shape::Tuple{Int64, Int64}, rng) = Dense(rand(rng, Float64, shape), zeros(shape[1]))

function layerEval(x, layer::Dense)
    return layer.W * x + layer.b
end


# --------------------- Activation Functions -----------------------
struct Exp <: Layer end
function layerEval(x, layer::Exp)
    return exp.(x)
end

struct DeArray <: Layer end
function layerEval(x, layer::DeArray)
    return sum(x)
end

struct Tanh <: Layer end
function layerEval(x, layer::Tanh)
    return tanh.(x)
end

struct Sigmoid <: Layer end
function layerEval(x, layer::Sigmoid)
    return 1 ./ (1 .+ exp.(-x))
end