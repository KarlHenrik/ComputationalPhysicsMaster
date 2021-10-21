abstract type Layer end

function layerParams(layer::Layer)
    return -1, -1
end

# --------------------- Dense -----------------------
struct Dense{M, V} <: Layer
    W::M
    b::V
end
Dense(shape::Tuple{Int64, Int64}, rng) = Dense(rand(rng, Float64, shape), zeros(shape[1]))

function layerEval(x, layer::Dense)
    return layer.W * x + layer.b
end

function layerEval(x, layer::Dense, params)
    W = params[1]
    b = params[2]
    return W * x + b
end

function layerParams(layer::Dense)
    flat = vcat(vec(layer.W), layer.b)
    shapes = [size(layer.W), size(layer.b)]
    return flat, shapes
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