# ------------Mutable storage for specific metropolis setups---------------

abstract type Metro_m end

struct Metro_m_basic <: Metro_m
    # Nothing needed
end

function get_metro_m(positions, wf, metro::Metropolis)
    return Metro_m_basic()
end

mutable struct Metro_m_imp <: Metro_m
    move::Vector{Float64}
    greens::Vector{Float64}
    
    new_pos::Vector{Float64}
    old_pos::Vector{Float64}
    oldQF::Vector{Float64}
    newQF::Vector{Float64}
end

function get_metro_m(positions, wf, metro::Importance)
    dims = system.dims
    return Metro_m_imp(zeros(dims), zeros(dims), zeros(dims), zeros(dims), zeros(dims), zeros(dims))
end

"""
mutable struct Metro_m_imp_slater <: Metro_m
    new_idx::UnitRange{Int64}
    new_pos::Vector{Float64}
    QF::Vector{Float64}
    ratio::Float64
    
    wf_der_mat::Matrix{Float64} # The derivative of each element in the slater determinant matrix
end

function get_metro_m(positions, system::System{W, H}, metro::Importance) where W <: Slater, H <: Hamiltonian
    return Metro_m_basic(0, 0, 0)
end
"""