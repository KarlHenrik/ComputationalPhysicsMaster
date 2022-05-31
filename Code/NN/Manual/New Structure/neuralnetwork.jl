struct NeuralNetwork{T}
    layers::Vector{Layer}
    output::Vector{Float64}
    
    # Parameter and input derivative, backpropogation
    delta_start::Vector{Float64}
    input_der::Vector{Float64}
    
    # Hessian
    hes_grad_result::Vector{Float64}
    hes_grad_tapes::Dict{DataType, Any}
    hes_jac_result::Matrix{Float64}
    hes_jac_config::T
end

function NeuralNetwork(layer_specs...; input, rng)
    layers = []
    
    output = zero(input)
    for spec in layer_specs
        layer_input = output # the input to the new layer is the output of the previous
        
        if typeof(spec) == DataType # Activation function
            layer, output = spec(layer_input)
        else # Dense layer
            spec, num_outputs = spec
            layer, output = spec(layer_input, num_outputs, rng)
        end
        
        push!(layers, layer)
    end
    
    # Parameter and input derivative setup
    delta_start = zero(layers[end].output)
    input_der = layers[1].delta
    
    # Kinetic energy (double derivative of network wrt. inputs)
    n = length(input)
    hes_grad_result = zeros(n)
    config = rd.GradientConfig(input)
    hes_grad_tapes = Dict{DataType, Any}()
    hes_jac_result = zeros(n, n)
    hes_jac_config = fd.JacobianConfig(nothing, hes_grad_result, hes_grad_result)
    
    return NeuralNetwork{typeof(hes_jac_config)}(
        layers, layers[end].output,
        delta_start, input_der,
        hes_grad_result, hes_grad_tapes, hes_jac_result, hes_jac_config
    )
end

function model!(nn::NeuralNetwork, x::Vector{Float64})
    (; layers) = nn
    layers[1].input .= x
    
    for layer in layers
        layerEval!(layer)
    end
    return layers[end].output
end

function gradient!(nn::NeuralNetwork)
    # Computes the derivative of the network output wrt. the parameters and input
    delta = nn.delta_start
    delta .= 1
    
    for layer in reverse(nn.layers)
        delta = backprop!(layer, delta)
    end
    
    return nn
end

function QF(nn)
    # TODO make this inplace?
    return 2 .* nn.input_der ./ nn.output[1]
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

function kinetic(x, nn::NeuralNetwork)
    # Assumes that model!(nn, x) has been called with the same x and nn as this call, since it uses the output of that computation
    fd.jacobian!(nn.hes_jac_result, (y, x) -> hes_grad!(y, x, nn), nn.hes_grad_result, x, nn.hes_jac_config)
    return -0.5 * la.tr(nn.hes_jac_result) / nn.output[1]
end
;