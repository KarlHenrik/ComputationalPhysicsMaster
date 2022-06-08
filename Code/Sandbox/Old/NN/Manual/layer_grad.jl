abstract type Layer_Grad end
setInput!(layer_grad, input) = return
setOutput!(layer_grad, output) = return
backprop!(layer_grad, delta, layer) = backprop!(layer_grad, delta)


# ------------------------ Dense ------------------------------
struct Dense_Grad <: Layer_Grad
    W_g::Matrix{Float64}
    b_g::Vector{Float64}
    
    input::Vector{Float64}
end
function Layer_Grad(layer::Dense, input)
    W_g = zero(layer.W)
    b_g = zero(layer.b)

    input = zero(input)

    return Dense_Grad(W_g, b_g, input)
end
function setInput!(layer_grad::Dense_Grad, input)
    layer_grad.input .= input
end
function backprop!(layer_grad::Dense_Grad, delta, layer)
    (; W_g, b_g, input) = layer_grad
    
    W_g .= delta * transpose(input)
    b_g .= delta
    
    delta = transpose(layer.W) * delta
    return delta
end


# ------------------------ Sigmoid ------------------------------
struct Sigmoid_Grad <: Layer_Grad
    output::Vector{Float64}
end
function Layer_Grad(layer::Sigmoid, input)
    return Sigmoid_Grad(zero(input))
end
function setOutput!(layer_grad::Sigmoid_Grad, output)
    layer_grad.output .= output
end
function backprop!(layer_grad::Sigmoid_Grad, delta)
    (; output) = layer_grad
    delta .= delta .* output .* (1 .- output)
end


# ------------------------ DeArray ------------------------------
struct DeArray_Grad <: Layer_Grad
    theOnes::Vector{Float64}
end
function Layer_Grad(layer::DeArray, input)
    return DeArray_Grad(zero(input) .+ 1.0)
end
function backprop!(layer_grad::DeArray_Grad, delta)
    delta = layer_grad.theOnes * delta
end