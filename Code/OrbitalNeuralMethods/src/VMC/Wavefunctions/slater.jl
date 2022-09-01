struct Slater <: WaveFunction
    dims::Int64
    num::Int64

    # This is the basis and C transform that defines how we evaluate the slater determinant
    C::Matrix{Float64} # lxl matrix from RHF which we use the num first columns of to produce our determinant
    basis::HObasis # The basis we will spin double

    # Here we store the candidate values for the matrices below
    new_eval::Vector{Float64} # new row in slater determinant matrix
    new_der::Vector{Float64}
    new_kin::Vector{Float64}

    # These hold the current state of the wavefunction. When a particle move is accepted, we update them
    amp_up::Fast_Det # n/2 x n/2 matrix that will be used to compute determinants when a spin-up particle is moved
    amp_down::Fast_Det # same as above, but for spin-down
    der_mat::Matrix{Float64} # n x n/2 matrix which holds the derivative of our slater determinant functions
    kin_mat::Matrix{Float64} # same as above, but double derivatives

    function Slater(num, dims, C, basis)
        @assert num%2 == 0
        @assert dims == 1

        amp_up, amp_down = setupSlaterMatrix(basis.base, C, num)
        

        return new(dims, num, C, basis, new_eval, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat)
    end
end

function Slater(state::RHFState)
    (; C, system) = state
    (; basis) = system
    return Slater(C, basis)
end

function Slater(system::System)
    (; basis, transform) = system
    C = transform
    @assert checkResticted(C)
    return Slater(C, basis)
end

function setupSlaterMatrix(basis::SpinBasis, C, positions)
    l = basis.base.l
    n = length(positions)

    #TODO these sizes are wrong!
    amp_up = zeros(n, n)
    amp_down = zeros(n, n)
    der_mat = zeros(n, n)
    kin_mat = zeros(n, n)

    for ϕ_i in 1:num # loop over the orbitals in the determinant
        for p_i in 1:num÷2
            amp = 0.0
            der = 0.0
            dder = 0.0
            for bs_i in 1:l # loop over the basis used to construct the orbitals 
                amp += C[bs_i, ϕ_i] * bs_eval[bs_i]
                der += C[bs_i, ϕ_i] * bs_der[bs_i]
                dder += C[bs_i, ϕ_i] * bs_dder[bs_i]
            end
            amp_up[ϕ, p_i] = 1
            amp_down[ϕ, p_i] = 1
            der_mat[ϕ, p_i] = 1
            kin_mat[ϕ, p_i] = 1
        end
    end
end

function computeNewRows!(wf, walker::Walker{S, M}, x::Float64)  where S <: Union{E_Muts, Grad_Muts} where M <: Imp_Muts
    (; basis, new_eval, new_der, new_dder) = wf

    bs_eval, bs_der, bs_dder = fast_ho_all!(x, basis)

    for ϕ_i in 1:num # loop over the orbitals in the determinant
        amp = 0.0
        der = 0.0
        dder = 0.0
        for bs_i in 1:l # loop over the basis used to construct the orbitals 
            amp += C[bs_i, ϕ_i] * bs_eval[bs_i]
            der += C[bs_i, ϕ_i] * bs_der[bs_i]
            dder += C[bs_i, ϕ_i] * bs_dder[bs_i]
        end
        new_eval[ϕ_i] = amp
        new_der[ϕ_i] = der
        new_dder[ϕ_i] = dder
    end
end

function setNewRows!(wf, walker::Walker{S, M}, new_idx::Int64)  where S <: Union{E_Muts, Grad_Muts} where M <: Metro_Muts
    (; num, new_eval, new_der, new_dder, amp_up, amp_down, der_mat, kin_mat) = wf
    if new_idx < num÷2
        change_row!(amp_up, new_eval, new_idx)
    else
        change_row!(amp_down, new_eval, new_idx-num÷2)
    end
    der_mat[new_idx, :] .= new_der
    kin_mat[new_idx, :] .= new_dder
end

function kinetic(positions, wf::Slater)::Float64
    (; amp_up, amp_down, kin_mat, num) = wf
    kin = 0.0
    for row in 1:num÷2
        kin += ratio_new_old_det(amp_up, kin_mat[row, :], row)
    end
    for (i, row) in enumerate(num÷2:num)
        kin += ratio_new_old_det(amp_down, kin_mat[row, :], i)
    end
    return kin
end

# Only used for the old qf
function QF!(qf, positions, idx, wf::SimpleGaussian)
    (; amp_up, amp_down, num) = wf

    if new_idx < num÷2
        qf .+= ratio_new_old_det(amp_up, der_mat[idx, :], idx)
    else
        qf .+= ratio_new_old_det(amp_down, der_mat[idx, :], idx-num÷2)
    end
    
    return qf
end

