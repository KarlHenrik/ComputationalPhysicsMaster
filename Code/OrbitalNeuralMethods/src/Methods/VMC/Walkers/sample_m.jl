# --------------Mutable storage for specific sampling setups-------------------

abstract type Sample_m end

struct Sample_m_pass <: Sample_m
    # Nothing is needed!
end

function get_sample_m(positions, wf, scheme::OneBody)
    return get_no_sample()
end

function get_no_sample()
    return Sample_m_pass()
end

mutable struct Gradient_m{T} <: Sample_m
    kinetic::Float64
    paramDer::T
end

function get_sample_m(positions, wf, scheme::OptimizerScheme)
    kinetic_ = kinetic(positions, wf)
    paramDer_ = paramDer(positions, wf)
    return Gradient_m(kinetic_, paramDer_)
end

mutable struct Blocking_m <: Sample_m
    kinetic::Float64
end

function get_sample_m(positions, wf, scheme::Blocking)
    kinetic_ = kinetic(positions, wf)
    return Blocking_m(kinetic_)
end

abstract type Sample_m_energy_slater <: Sample_m end

"""
mutable struct Sample_m_gradient_slater{T} <: Sample_m_energy_slater
    kinetic::Float64
    paramDer::T
    wf_dder_mat::Matrix{Float64}
    wf_dder_newcol::Vector{Float64}
end

mutable struct Sample_m_blocking_slater <: Sample_m_energy_slater
    kinetic::Float64
    wf_dder_mat::Matrix{Float64}
    wf_dder_newcol::Vector{Float64}
end
"""