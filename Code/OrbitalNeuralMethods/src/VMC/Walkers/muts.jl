struct No_Muts end

mutable struct E_Muts
    kinetic::Float64
    potential::Float64
    E::Float64
    E2::Float64
end

mutable struct Grad_Muts{T}
    kinetic::Float64
    potential::Float64
    E::Float64
    E2::Float64
    paramDer::T
end
Sample_Muts(wf, sampler::OneBodySampler) = No_Muts()
Sample_Muts(wf, sampler::BlockingSampler) = E_Muts(0, 0, 0, 0)
function Sample_Muts(wf, sampler::GradientSampler)
    paramDer = paramDerHolder(wf)
    return Grad_Muts{typeof(paramDer)}(0, 0, 0, 0, paramDer)
end