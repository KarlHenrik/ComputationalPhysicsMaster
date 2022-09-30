struct NeuralNetwork{N}
    n::Int64
    layers::Vector{Layer}
    output::Vector{Float64}
    old_output::Vector{Float64}
    
    # Parameter and input derivative, backpropogation
    delta_start::Vector{Float64}
    input_der::Vector{Float64}
    QF_all_old::Vector{Float64}
    
    # Hessian
    hes_grad_result::Vector{Float64}
    hes_grad_tapes::Dict{DataType, Any}
    hes_jac_result::Matrix{Float64}
    hes_jac_config::fd.JacobianConfig{Nothing, Float64, N, Tuple{Vector{fd.Dual{Nothing, Float64, N}}, Vector{fd.Dual{Nothing, Float64, N}}}}
end

function NeuralNetwork(layers::Vector{Layer}, n::Int64)
    # Parameter and input derivative setup
    delta_start = zero(layers[end].output)
    input_der = layers[1].delta
    QF_all_old = zero(input_der)
    
    # Kinetic energy (double derivative of network wrt. inputs)
    hes_grad_result = zeros(n)
    hes_grad_tapes = Dict{DataType, Any}()
    hes_jac_result = zeros(n, n)
    hes_jac_config = fd.JacobianConfig(nothing, hes_grad_result, hes_grad_result)
    
    return NeuralNetwork{chunksize(fd.Chunk(hes_grad_result))}(
        n, layers, layers[end].output, zero(layers[end].output),
        delta_start, input_der, QF_all_old,
        hes_grad_result, hes_grad_tapes, hes_jac_result, hes_jac_config
    )
end

function chunksize(chunk::fd.Chunk{N}) where N
    return N
end

function NeuralNetwork(layer_specs; n::Int64, rng)
    layers = []

    output = zeros(n)
    for spec in layer_specs
        layer_input = output # the input to the new layer is the output of the previous
        
        if typeof(spec) == DataType # Activation function
            layer, output = spec(layer_input)
        else # Dense
            spec, num_outputs = spec
            layer, output = spec(layer_input, num_outputs, rng)
        end
        
        push!(layers, layer)
    end
    layers = Vector{Layer}(layers)
    return NeuralNetwork(layers, n)
end

function private_wf(nn::NeuralNetwork, positions)
    n = nn.n
    layers = []
    output = zeros(n)
    for layer in nn.layers
        layer_input = output # the input to the new layer is the output of the previous
        
        layer, output = layerCopy(layer_input, layer)
        
        push!(layers, layer)
    end
    layers = Vector{Layer}(layers)
    
    wf = NeuralNetwork(layers, n)
    model!(wf, positions)
    gradient!(wf)
    wf.old_output[1] = wf.output[1]
    wf.QF_all_old .= 2 .* wf.input_der ./ wf.output[1]
    return wf
end

function model!(nn::NeuralNetwork, x::Vector{Float64})::Vector{Float64}
    (; layers) = nn
    layers[1].input::Vector{Float64} .= x
    
    for layer in layers
        layerEval!(layer)
    end
    return layers[end].output
end

function gradient!(nn::NeuralNetwork{T}) where T
    # Computes the derivative of the network output wrt. the parameters and input
    delta = nn.delta_start
    delta .= 1
    
    for layer in Iterators.reverse(nn.layers)
        delta = backprop!(layer, delta)::Vector{Float64}
    end
    
    return nn::NeuralNetwork{T}
end

struct Dense_Grad
    layer_idx::Int64
    W_g::Matrix{Float64}
    b_g::Vector{Float64}
end
function Layer_Grad(layer::Dense, i)
    return Dense_Grad(i, zero(layer.W_g), zero(layer.b_g))
end

function getGrad!(layer_grad::Dense_Grad, layer::Dense)
    layer_grad.W_g .= layer.W_g
    layer_grad.b_g .= layer.b_g

    return layer_grad
end

function add!(d1::Dense_Grad, d2::Dense_Grad)
    d1.W_g .+= d2.W_g
    d1.b_g .+= d2.b_g
    return d1
end
function add!(grads1::Vector{Dense_Grad}, grads2::Vector{Dense_Grad})
    for (d1, d2) in zip(grads1, grads2)
        add!(d1, d2)
    end
    return grads1
end

import Base.*
function *(d::Dense_Grad, f::Number)
    d.W_g .*= f
    d.b_g .*= f
    return d
end
import Base./
function /(d::Dense_Grad, f::Number)
    d.W_g ./= f
    d.b_g ./= f
    return d
end
function setmul!(d1::Dense_Grad, d2::Dense_Grad, f::Float64)
    d1.W_g .= d2.W_g .* f
    d1.b_g .= d2.b_g .* f

    return d1
end
function setmul!(grads1::Vector{Dense_Grad}, grads2::Vector{Dense_Grad}, f::Float64)
    for (d1, d2) in zip(grads1, grads2)
        setmul!(d1, d2, f)
    end
    return grads1
end

function paramDerHolder(nn::NeuralNetwork)
    layer_grads = []
    for (i, layer) in enumerate(nn.layers)
        if typeof(layer) == Dense
            push!(layer_grads, Layer_Grad(layer, i))
        end
    end
    layer_grads = Vector{Dense_Grad}(layer_grads)
    return layer_grads
end
function la.norm(layer_grads::Vector{Dense_Grad})
    sq_sum = 0.0
    for layer_grad in layer_grads
        sq_sum += sum(layer_grad.W_g.^2)
        sq_sum += sum(layer_grad.b_g.^2)
    end
    return √sq_sum
end

function paramDer!(layer_grads::Vector{Dense_Grad}, positions, wf::NeuralNetwork)
    amp = wf.output[1]
    layers = wf.layers
    for layer_grad in layer_grads
        layer_grad = getGrad!(layer_grad, layers[layer_grad.layer_idx]::Dense)
        layer_grad = layer_grad / amp
    end
    return layer_grads
end

function applyGradient(nn::NeuralNetwork, grad::Vector{Dense_Grad})
    n = nn.n
    layers = []
    output = zeros(n)
    i = 1
    for layer in nn.layers
        layer_input = output # the input to the new layer is the output of the previous
        
        layer, output = layerCopy(layer_input, layer)
        if typeof(layer) == Dense
            layer_grad = grad[i]
            i += 1
            layer.W .-= layer_grad.W_g
            layer.b .-= layer_grad.b_g
        end
        
        push!(layers, layer)
    end
    layers = Vector{Layer}(layers)
    return NeuralNetwork(layers, n)
end

# Only used for old QF
function QF(positions, new_idx, wf)
    return wf.QF_all_old[new_idx]
end

function diffable_model(layers, x)
    for layer in layers
        x = layerEval(x, layer)
    end
    return x
end

function hes_grad!(y, x::Array{T}, nn::NeuralNetwork) where {T<:Real}
    if !haskey(nn.hes_grad_tapes, T)
        config = rd.GradientConfig(x)
        tape = rd.compile(rd.GradientTape(x -> diffable_model(nn.layers, x)[1], x, config))
        nn.hes_grad_tapes[T] = tape
    end
    tape = nn.hes_grad_tapes[T]
    return rd.gradient!(y, tape, x)
end

function kinetic(x::Vector{Float64}, nn::NeuralNetwork)
    # Assumes that model!(nn, x) has been called with the same x and nn as this call, since it uses the output of that computation
    fd.jacobian!(nn.hes_jac_result, (y, x) -> hes_grad!(y, x, nn), nn.hes_grad_result, x, nn.hes_jac_config)
    return -0.5 * la.tr(nn.hes_jac_result) / nn.old_output[1]
end

function dder!(nn_dder, x, nn::NeuralNetwork)
    # Assumes that model!(nn, x) has been called with the same x and nn as this call, since it uses the output of that computation
    fd.jacobian!(nn.hes_jac_result, (y, x) -> hes_grad!(y, x, nn), nn.hes_grad_result, x, nn.hes_jac_config)
    for i in eachindex(nn_dder)
        nn_dder[i] = nn.hes_jac_result[i, i]
    end
    nn_dder .= nn_dder ./ nn.old_output[1]
    return nn_dder
end

# Importance + NeuralNetwork + ParamDer is set up
function consider_qf!(wf::NeuralNetwork, positions, new_idx::Int64, old_pos)
    output = model!(wf, positions)
    gradient!(wf)
    newQF = 2 * wf.input_der[new_idx] / output[1]
    
    ratio = output[1] / wf.old_output[1]
    return ratio, newQF
end

# Importance + NeuralNetwork + ParamDer is set up
function consider!(wf::NeuralNetwork, positions, new_idx::Int64, old_pos)
    output = model!(wf, positions)
    gradient!(wf)
    
    ratio = output[1] / wf.old_output[1]
    return ratio
end

function accept!(wf::NeuralNetwork, new_idx, new_pos)
    wf.old_output[1] = wf.output[1]
    wf.QF_all_old .= 2 .* wf.input_der ./ wf.output[1]
    
    return wf
end
;