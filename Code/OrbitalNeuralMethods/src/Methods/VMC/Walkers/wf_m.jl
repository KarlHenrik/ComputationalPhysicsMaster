# ---------------Mutable storage for specific wavefunctions-------------------

abstract type Wf_m end

struct Wf_m_pass <: Wf_m
    # Nothing is needed! 
end

function get_wf_m(positions, wf::WaveFunction)
    return Wf_m_pass()
end

"""

mutable struct Wf_m_slater
    amplitude::Float64
    wf_mat::Matrix{Float64} # The elements in the slater determinant matrix
end

function get_wf_m(positions, system{W, H}) where W <: Slater, H <: Hamiltonian
    # TODO
    wf_mat = 
    amplitude = 
    return Wf_m_pass(amplitude, wf_mat)
end
"""