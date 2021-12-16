struct System{W<:WaveFunction, H<:Hamiltonian}
    dims::Int64
    num::Int64
    wf::W
    ham::H
end