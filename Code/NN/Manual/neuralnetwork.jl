

struct NeuralNetwork
    layers::Vector{Layer}
    layer_grads::Vector{Layer_Grad}
    
    function NeuralNetwork(layers::Vector{Layer}, x)
        input = x
        layer_grads = Vector{Layer_Grad}(undef, length(layers))
        
        for (i, layer) in enumerate(layers)
            layer_grads[i] = Layer_Grad(layer, input)
            input = layerEval(input, layer)
        end
        
        return new(layers, layer_grads)
    end
end

function model(x, nn::NeuralNetwork)
    for layer in nn.layers
        x = layerEval(x, layer)
    end
    return x
end

function kinetic!(nn::NeuralNetwork, x)
    return
end

function gradient!(nn::NeuralNetwork, x)
    for (layer, layer_grad) in zip(nn.layers, nn.layer_grads)
        setInput!(layer_grad, x)
        x = layerEval(x, layer)
        setOutput!(layer_grad, x)
    end
    
    delta = x
    for (layer, layer_grad) in zip(reverse(nn.layers), reverse(nn.layer_grads))
        backprop!(layer_grad, delta, layer)
    end
    
    return nn
end