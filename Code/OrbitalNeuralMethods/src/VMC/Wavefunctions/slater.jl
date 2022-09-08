struct Slater <: WaveFunction
    dims::Int64
    num::Int64

    # This is the basis and C transform that defines how we evaluate the slater determinant
    C::Matrix{Float64} # lxl matrix from RHF which we use the num first columns of to produce our determinant
    basis::HOBasis # The basis we will spin double

    # Here we store the candidate values for the matrices below
    new_amp::Vector{Float64} # new row in slater determinant matrix
    new_der::Vector{Float64}
    new_kin::Vector{Float64}

    # These hold the current state of the wavefunction. When a particle move is accepted, we update them
    amp_up::Fast_Det # n/2 x n/2 matrix that will be used to compute determinants when a spin-up particle is moved
    amp_down::Fast_Det # same as above, but for spin-down
    der_mat::Matrix{Float64} # n x n/2 matrix which holds the derivative of our slater determinant functions
    kin_mat::Matrix{Float64} # same as above, but double derivatives
end

# * This is the "true" constructor, since the saved values depend on the positions of the electons 
function private_wf(wf::Slater, positions)
    (; dims, num, C, basis) = wf
    amp_up, amp_down, der_mat, kin_mat = setupSlaterMatrix(basis, C, positions)
    new_amp, new_der, new_kin = zeros(num÷2), zeros(num÷2), zeros(num÷2)
    #! Is this right?
    return Slater(dims, num, C, basis, new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat)
end

# * These other constructors only make sure the basis and C are correct
function Slater(num, C, basis::SpinBasis)
    @assert num%2 == 0
    dims = 1
    # These are placeholder values, private_wf gives an actual usable Slater
    amp_up = Fast_Det(Matrix(Float64.(la.I(num÷2))))
    amp_down = Fast_Det(Matrix(Float64.(la.I(num÷2))))
    new_amp, new_der, new_kin = zeros(num), zeros(num), zeros(num)
    der_mat, kin_mat = zeros(num, num÷2), zeros(num, num÷2)

    return Slater(dims, num, C, basis.base, new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat)
end

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

function Slater(num, basis::SpinBasis)
    (; l) = basis.base
    @assert num%2 == 0
    C = Matrix(Float64.(la.I(l)))

    return Slater(num, C, basis)
end

function setupSlaterMatrix(basis::HOBasis, C, positions)
    l = basis.l
    n = length(positions)

    amp_up = zeros(n÷2, n÷2)
    amp_down = zeros(n÷2, n÷2)
    der_mat = zeros(n, n÷2)
    kin_mat = zeros(n, n÷2)

    for r_i in 1:n # loop over the particles in the determiant
        bs_eval, bs_der, bs_dder = fast_ho_all!(positions[r_i], basis)
        for ϕ_i in 1:n÷2 # loop over the orbitals in the determinant
            amp = 0.0
            der = 0.0
            dder = 0.0
            for bs_i in 1:l # loop over the basis used to construct the orbitals 
                amp += C[bs_i, ϕ_i] * bs_eval[bs_i]
                der += C[bs_i, ϕ_i] * bs_der[bs_i]
                dder += C[bs_i, ϕ_i] * bs_dder[bs_i]
            end
            if r_i <= n÷2 # The spin-up particles
                amp_up[r_i, ϕ_i] = amp

                der_mat[r_i, ϕ_i] = der
                kin_mat[r_i, ϕ_i] = dder
            else         # The spin-down particles
                amp_down[r_i - n÷2, ϕ_i] = amp

                der_mat[r_i, ϕ_i] = der
                kin_mat[r_i, ϕ_i] = dder
            end
        end
    end
    fast_amp_up = Fast_Det(amp_up)
    fast_amp_down = Fast_Det(amp_down)
    return fast_amp_up, fast_amp_down, der_mat, kin_mat
end

function computeNewRows!(wf::Slater, x::Float64)
    (; basis, new_amp, new_der, new_kin, num, C) = wf
    l = size(C)[1]
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
        #! Figure out the shape of the derivatives and amplitues
        new_amp[ϕ_i] = amp
        new_der[ϕ_i] = der
        new_kin[ϕ_i] = dder
    end
    return wf
end

#* This only makes sence after computeNewRows! has been called
function ratio_direct(wf::Slater, new_idx)
    (; num, amp_up, amp_down, new_amp) = wf

    if new_idx <= num÷2
        ratio = ratio_new_old_det(amp_up, new_amp, new_idx)
    else
        ratio = ratio_new_old_det(amp_down, new_amp, new_idx-num÷2)
    end
    return ratio
end

function setNewRows!(wf::Slater, new_idx::Int64)
    (; num, new_amp, new_der, new_kin, amp_up, amp_down, der_mat, kin_mat) = wf
    if new_idx <= num÷2
        change_row!(amp_up, new_amp, new_idx)
    else
        change_row!(amp_down, new_amp, new_idx-num÷2)
    end
    println(new_der)
    println(der_mat)
    der_mat[new_idx, :] .= new_der
    kin_mat[new_idx, :] .= new_kin
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

#* Only used for the old qf
function QF!(qf, positions, idx, wf::Slater)
    (; amp_up, amp_down, der_mat, num) = wf
    new_idx = idx[1]
    if new_idx <= num÷2
        qf .+= ratio_new_old_det(amp_up, der_mat[new_idx, :], new_idx)
    else
        qf .+= ratio_new_old_det(amp_down, der_mat[new_idx, :], new_idx-num÷2)
    end
    
    return qf
end

function paramDer!(samp_muts, positions, wf::Slater)
    throw("Slater determinants should not be used for gradients")
end