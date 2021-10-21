# ----------------------- Utility for working with flat parameters ----------------
struct Unflattener
    indices::Vector{UnitRange{Int64}}
    shapes::Vector{Tuple{Int64, Vararg{Int64, N} where N}}
    num_params::Int
end

function unflatten(unflat, params_flat)
    return [reshape(params_flat[unflat.indices[i]], unflat.shapes[i]) for i in 1:unflat.num_params]
end

# -----------------------------------------------------------------------------------

function flattenNetwork(layers)
    params_flat = []
    param_unflatteners = []
    index = 1
    for layer in layers
        flat, shapes = layerParams(layer)
        indices = []
        if shapes != -1
            params_flat = vcat(params_flat, flat)
            for i in 1:length(shapes)
                push!(indices, index:index + prod(shapes[i]) - 1)
                index += prod(shapes[i])
            end
            num_params = length(shapes)
        else
            shapes = [(0,)]
            indices = [2:1]
            num_params = 0
        end
        push!(param_unflatteners, Unflattener(indices, shapes, num_params))
    end
    return Float64.(params_flat), param_unflatteners
end

struct NeuralNetwork{T, C}
    layers::Vector{Layer}
    # Flat representation, and unflattener. Used for taking the gradient of the network parameters
    params_flat::Vector{Float64}
    unflatteners::Vector{Unflattener}
    # Gradient setup
    grad_result::Vector{Float64}
    grad_tape::T
    # Hessian setup
    hes_grad_result::Vector{Float64}
    hes_grad_tapes::Dict{DataType, Any}
    hes_jac_result::Matrix{Float64}
    hes_jac_config::C
    
    function NeuralNetwork(layers::Vector{Layer}, x)
        params_flat, unflatteners = flattenNetwork(layers)
        
        n = length(x)
        grad_input = vcat(x, params_flat)
        grad_result = zero(grad_input)
        grad_config = rd.GradientConfig(grad_input)
        grad_tape = rd.GradientTape(grad_input -> model_flat(grad_input[1:n], grad_input[n+1:end], layers, unflatteners), grad_input, grad_config)
        grad_tape = rd.compile(grad_tape)
        
        hes_grad_result = zeros(n)
        config = rd.GradientConfig(x)
        hes_grad_tapes = Dict{DataType, Any}()
        hes_jac_result = zeros(n, n)
        hes_jac_config = fd.JacobianConfig(nothing, hes_grad_result, hes_grad_result)
        
        return new{typeof(grad_tape), typeof(hes_jac_config)}(layers, params_flat, unflatteners,
                                                              grad_result, grad_tape,
                                                              hes_grad_result, hes_grad_tapes, hes_jac_result, hes_jac_config)
    end
end

function model(x, nn::NeuralNetwork)
    for layer in nn.layers
        x = layerEval(x, layer)
    end
    return x
end

function model_flat(x, nn::NeuralNetwork, params_flat)
    for j in 1:length(nn.layers)
        unflat = nn.unflatteners[j]
        if unflat.num_params == 0
            x = layerEval(x, nn.layers[j])
        else
            x = layerEval(x, nn.layers[j], unflatten(unflat, params_flat))
        end
    end
    return x
end

function model_flat(x, params_flat, layers, unflatteners)
    for j in 1:length(layers)
        unflat = unflatteners[j]
        if unflat.num_params == 0
            x = layerEval(x, layers[j])
        else
            x = layerEval(x, layers[j], unflatten(unflat, params_flat))
        end
    end
    return x
end

function hes_grad!(y, x::Array{T}, nn::NeuralNetwork) where {T<:Real}
    if !haskey(nn.hes_grad_tapes, T)
        config = rd.GradientConfig(x)
        tape = rd.compile(rd.GradientTape(x -> model(x, nn), x, config))
        nn.hes_grad_tapes[T] = tape
    end
    tape = nn.hes_grad_tapes[T]
    return rd.gradient!(y, tape, x)
end

function kinetic(x, nn::NeuralNetwork)
    fd.jacobian!(nn.hes_jac_result, (y, x) -> hes_grad!(y, x, nn), nn.hes_grad_result, x, nn.hes_jac_config)
    return -0.5 * la.tr(nn.hes_jac_result)
end

function gradient(x, nn::NeuralNetwork)
    result = rd.gradient!(nn.grad_result, nn.grad_tape, vcat(x, nn.params_flat))
    return result[length(x)+1:end]
end