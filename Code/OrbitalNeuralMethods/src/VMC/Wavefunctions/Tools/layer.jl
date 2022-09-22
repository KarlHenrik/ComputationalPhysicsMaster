abstract type Layer end


# --------------------- Dense -----------------------
struct Dense <: Layer
    W::Matrix{Float64}
    b::Vector{Float64}
    
    input::Vector{Float64}
    output::Vector{Float64}
    
    W_g::Matrix{Float64}
    b_g::Vector{Float64}
    delta::Vector{Float64}
end

function Dense(input::Vector{Float64}, W::Matrix{Float64}, b::Vector{Float64})
    output = zero(b)
    
    W_g = zero(W)
    b_g = zero(b)
    delta = zero(input)
    
    return Dense(W, b, input, output, W_g, b_g, delta), output
end

# Initial construction
function Dense(input::Vector{Float64}, num_outputs::Int64, rng)
    n = length(input)
    
    W = randn(rng, Float64, (num_outputs, n)) ./ n
    #b = randn(rng, Float64, num_outputs) ./ n
    b = zeros(num_outputs)
    
    Dense(input, W, b)
end

function layerCopy(input::Vector{Float64}, layer::Dense)
    W = copy(layer.W)
    b = copy(layer.b)
    
    Dense(input, W, b)
end

Dense(num_outputs::Int) = (Dense, num_outputs)

function layerEval!(layer::Dense)
    (; W, b, input, output) = layer
    
    la.mul!(output, W, input)
    output .+= b
end

function layerEval(x, layer::Dense)
    (; W, b) = layer
    
    return W * x + b
end

#* Faster when the dimentions are smaller than about 30
function vec_inner!(M::Matrix{Float64}, u, vT)
    @inbounds for (i, ui) in enumerate(u)
        for (j, vi) in enumerate(vT)
            M[i, j] = ui * vi
        end
    end
    return M
end

function backprop!(layer::Dense, delta_in::Vector{Float64})
    (; W, W_g, b_g, input, delta) = layer
    
    #la.mul!(W_g, delta_in, transpose(input)) # This is faster if we have many nodes in each layer
    vec_inner!(W_g, delta_in, input) # Weight gradient
    b_g .= delta_in # Bias gradient
    
    la.mul!(delta, transpose(W), delta_in) # Propagating gradient
    return delta
end


# --------------------- Activation Functions -----------------------
abstract type Activation <: Layer end

function layerCopy(input::Vector{Float64}, layer::Activation)
    return typeof(layer)(input)
end

#Exp
struct Exp <: Activation
    input::Vector{Float64}
    output::Vector{Float64}
end
function Exp(input)
    output = zero(input)
    return Exp(input, output), output
end
function layerEval(x, layer::Exp)
    return exp.(x)
end
function layerEval!(layer::Exp)
    (; input, output) = layer
    output .= exp.(input)
end
function backprop!(layer::Exp, delta)
    (; output) = layer
    delta .*= output
end


# Tanh
struct Tanh <: Activation
    input::Vector{Float64}
    output::Vector{Float64}
end
function Tanh(input)
    output = zero(input)
    return Tanh(input, output), output
end
function layerEval(x, layer::Tanh)
    return tanh.(x)
end
function layerEval!(layer::Tanh)
    (; input, output) = layer
    output .= tanh.(input)
end
function backprop!(layer::Tanh, delta)
    (; input) = layer
    delta .*= sech.(input).^2
end




# Sigmoid
function sig(x)
    return 1.0 ./ (1.0 .+ exp.(-x)) # TODO do i need two definitions? with a "where"?
end
struct Sigmoid <: Activation
    input::Vector{Float64}
    output::Vector{Float64}
end
function Sigmoid(input)
    output = zero(input)
    return Sigmoid(input, output), output
end
function layerEval(x, layer::Sigmoid)
    return sig.(x)
end
function layerEval!(layer::Sigmoid)
    (; input, output) = layer
    
    output .= sig.(input)
end
function backprop!(layer::Sigmoid, delta)
    (; output) = layer
    delta .*= output .* (1 .- output)
end
;