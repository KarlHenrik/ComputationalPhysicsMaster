struct System{W<:WaveFunction, H<:Hamiltonian}
    dims::Int64
    num::Int64
    wf::W
    ham::H
    function System(wf, ham)
        @assert wf.dims == ham.dims
        @assert wf.num  == ham.num
        return new{typeof(wf), typeof(ham)}(wf.dims, wf.num, wf, ham)
    end
end