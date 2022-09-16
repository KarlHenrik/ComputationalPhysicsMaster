struct Slater <: WaveFunction
    n::Int64

    # This is the basis and C transform that defines how we evaluate the slater determinant
    C::Matrix{Float64} # l x n/2 matrix from RHF which we use the columns of to produce our determinant
    basis::HOBasis # The basis we will spin double

    # In-place evaluation of the basis
    basis_tmp::Vector{Float64}
    basis_amp::Vector{Float64}
    basis_der::Vector{Float64}
    basis_kin::Vector{Float64}

    # Here we store the candidate values for the matrices below
    new_amp::Vector{Float64} # new row in slater determinant matrix
    new_der::Vector{Float64} # these have length n/2
    new_kin::Vector{Float64}

    # These hold the current state of the wavefunction. When a particle move is accepted, we update them
    amp_up::Fast_Det # n/2 x n/2 matrix that will be used to compute determinants when a spin-up particle is moved
    amp_down::Fast_Det # same as above, but for spin-down
    der_mat::Matrix{Float64} # n x n/2 matrix which holds the derivative of our slater determinant functions
    kin_mat::Matrix{Float64} # same as above, but double derivatives
end

#* INITIALIZES A USABLE SLATER, since the saved values depend on the positions of the electons 
function initSlater(n, C, basis::HOBasis, positions)
    basis_tmp, basis_amp, basis_der, basis_kin = zeros(basis.l), zeros(basis.l), zeros(basis.l), zeros(basis.l)

    amp_up, amp_down, der_mat, kin_mat = setupSlaterMatrix(basis, C, positions, basis_tmp, basis_amp, basis_der, basis_kin)

    new_amp, new_der, new_kin = zeros(n÷2), zeros(n÷2), zeros(n÷2)
    
    return Slater(n, C, basis, basis_tmp, basis_amp, basis_der, basis_kin,
                    new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat)
end

function private_wf(wf::Slater, positions)
    (; n, C, basis) = wf
    return initSlater(n, C, basis, positions)
end

#* MAIN CONSTRUCTOR OF AN UNUSABLE SLATER TEMPLATE: Sets up n and the required C columns
# C must be only lxl, where l is the number of unique spatial functions
# The spin basis requirement is for theoretical coherence, we actually include spin, even if we don't use the spinbasis
function Slater(n, C, basis::SpinBasis)
    @assert n%2 == 0

    C = C[:, 1:n÷2]

    # These are placeholder values, private_wf gives an actual usable Slater
    basis_tmp, basis_amp, basis_der, basis_kin = zeros(0), zeros(0), zeros(0), zeros(0)
    amp_up, amp_down, der_mat, kin_mat = Fast_Det(zeros(0,0)), Fast_Det(zeros(0,0)), zeros(0, 0), zeros(0, 0)
    new_amp, new_der, new_kin = zeros(0), zeros(0), zeros(0)

    return Slater(n, C, basis.base, basis_tmp, basis_amp, basis_der, basis_kin,
                    new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat)
end

# Different ways of making a Slater template
function Slater(state::RHFState)
    (; C, system) = state
    (; basis, n) = system
    @assert n%2 == 0
    return Slater(n, C, basis)
end

function Slater(system::SpatialSystem{SpinBasis})
    (; basis, transform, n) = system
    @assert n%2 == 0

    C = transform
    @assert checkResticted(C)
    C = shrink_restricted(C)

    return Slater(n, C, basis)
end

function Slater(n, basis::SpinBasis)
    (; l) = basis.base
    @assert n%2 == 0
    C = Matrix(Float64.(la.I(l)))

    return Slater(n, C, basis)
end



function setupSlaterMatrix(basis::HOBasis, C, positions, basis_tmp, basis_amp, basis_der, basis_kin)
    l = basis.l
    n = length(positions)

    amp_up = zeros(n÷2, n÷2)
    amp_down = zeros(n÷2, n÷2)
    der_mat = zeros(n, n÷2)
    kin_mat = zeros(n, n÷2)

    for r_i in 1:n # loop over the particles in the determiant
        basis_amp, basis_der, basis_kin = fast_ho_all!(basis_tmp, basis_amp, basis_der, basis_kin, positions[r_i], basis)
        for ϕ_i in 1:n÷2 # loop over the orbitals in the determinant
            amp = 0.0
            der = 0.0
            kin = 0.0
            for bs_i in 1:l # loop over the basis used to construct the orbitals 
                amp += C[bs_i, ϕ_i] * basis_amp[bs_i]
                der += C[bs_i, ϕ_i] * basis_der[bs_i]
                kin += C[bs_i, ϕ_i] * basis_kin[bs_i]
            end
            if r_i <= n÷2 # The spin-up particles
                amp_up[r_i, ϕ_i] = amp
            else         # The spin-down particles
                amp_down[r_i - n÷2, ϕ_i] = amp
            end
            der_mat[r_i, ϕ_i] = der
            kin_mat[r_i, ϕ_i] = kin
        end
    end
    fast_amp_up = Fast_Det(amp_up)
    fast_amp_down = Fast_Det(amp_down)
    return fast_amp_up, fast_amp_down, der_mat, kin_mat
end

function computeNewRows!(wf::Slater, x::Float64)
    (; basis, basis_tmp, basis_amp, basis_der, basis_kin, new_amp, new_der, new_kin, n, C) = wf
    l = size(C)[1]
    basis_amp, basis_der, basis_kin = fast_ho_all!(basis_tmp, basis_amp, basis_der, basis_kin, x, basis)

    for ϕ_i in 1:n÷2 # loop over the orbitals in the determinant
        amp = 0.0
        der = 0.0
        kin = 0.0
        for bs_i in 1:l # loop over the basis used to construct the orbitals 
            amp += C[bs_i, ϕ_i] * basis_amp[bs_i]
            der += C[bs_i, ϕ_i] * basis_der[bs_i]
            kin += C[bs_i, ϕ_i] * basis_kin[bs_i]
        end
        new_amp[ϕ_i] = amp
        new_der[ϕ_i] = der
        new_kin[ϕ_i] = kin
    end
    return wf
end

#* This only makes sence after computeNewRows! has been called
function ratio_direct(wf::Slater, new_idx)
    (; n, amp_up, amp_down, new_amp) = wf

    if new_idx <= n÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-n÷2)
    end
    return ratio
end

function setNewRows!(wf::Slater, new_idx::Int64)
    (; n, new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat) = wf
    if new_idx <= n÷2
        change_row!(amp_up, new_amp, new_idx)
    else
        change_row!(amp_down, new_amp, new_idx-n÷2)
    end
    der_mat[new_idx, :] .= new_der
    kin_mat[new_idx, :] .= new_kin

    return wf
end

function kinetic(positions, wf::Slater)::Float64
    (; amp_up, amp_down, kin_mat, n) = wf
    kin = 0.0
    #TODO maybe find faster way of doing this than slices
    for row in 1:n÷2
        @views kin += ratio_new_old_det(amp_up, kin_mat[row, :], row)
    end
    for (i, row) in enumerate(n÷2+1:n)
        @views kin += ratio_new_old_det(amp_down, kin_mat[row, :], i)
    end
    return -0.5 * kin
end

#* Only used for the old qf
function QF(positions, new_idx::Int64, wf::Slater)
    (; amp_up, amp_down, der_mat, n) = wf

    if new_idx <= n÷2
        @views qf = ratio_new_old_det(amp_up, der_mat[new_idx, :], new_idx)
    else
        @views qf = ratio_new_old_det(amp_down, der_mat[new_idx, :], new_idx-n÷2)
    end
    
    return qf
end

function paramDer!(samp_muts, positions, wf::Slater)
    throw("Slater determinants should not be used for gradients")
end

# ------------------- No first or second derivative ------------------------
function consider!(wf::Slater, positions, new_idx::Int64, old_pos)
    wf = computeNewRows!(wf, positions[new_idx])
    (; amp_up, amp_down, new_amp, n) = wf
    
    if new_idx <= n÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-n÷2)
    end

    return ratio
end

function consider_qf!(wf::Slater, positions, new_idx::Int64, old_pos)
    wf = computeNewRows!(wf, positions[new_idx])
    (; amp_up, amp_down, new_amp, new_der, n) = wf
    
    if new_idx <= n÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
        newQF = ratio_new_old_det(amp_up, new_der, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-n÷2)
        newQF = ratio_new_old_det(amp_down, new_der, new_idx-n÷2)
    end

    return ratio, newQF
end

function accept!(wf::Slater, new_idx::Int64)
    wf = setNewRows!(wf, new_idx)
    return wf
end