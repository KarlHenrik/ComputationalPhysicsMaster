# QUANTUM DYNAMICS
# by Jonas Boym Flaten



@enum ParticleClass boson fermion # encodes the two classes of particles.

struct Particle
    # represents a particle through its three main properties.

    mass::Float64       # is the particle's mass in mₑ.
    charge::Int64       # is the particle's charge in qₑ.
    spin::Rational      # is the particle's spin number in ħ.
        class::ParticleClass  # is the particle's class, inferred from 'spin'.

    function Particle(mass,charge,spin)

        @assert (mass ≥ 0.0) "The given mass $mass is not valid; input 0 or a positive number as 'mass'."

        @assert (spin%(1//2) == 0) "The given spin number $spin is not valid; input 0 or a positive multiple of 1/2 as 'spin'."

        if (spin%1 == 0)
            class = boson
        else
            class = fermion
        end

        return new(mass,charge,spin,class)
    end
end

Electron() = Particle(1.0,-1,1//2) # returns a Particle encoding an electron.
Proton() = Particle(1836.15,+1,1//2) # returns a Particle encoding a proton.
Neutron() = Particle(1838.68,0,1//2) # returns a Particle encoding a neutron.

Electrons(n) = (n,Electron()) # returns a tuple encoding n electrons.
Protons(n) = (n,Proton()) # returns a tuple encoding n protons.
Neutrons(n) = (n,Neutron()) # returns a tuple encoding n neutrons.



include("EnergyTerms.jl") # includes the various possible energy terms of the system Hamiltonian defined in EnergyTerms.jl.



struct System
    # represents a system through its content and dimensions.

    particles::Vector{Tuple{Int,Particle}} # are the particle groups of the system, given in tuples of particle number and type.
        N::Int64 # is the total number of particles in the system, inferred from 'particles'.
    Hamiltonian::Vector{EnergyTerm} # is the Hamiltonian of the system, given as a vector of EnergyTerms.
    dimensions::Vector{NTuple{2,Float64}} # are the dimensions of the system, given in tuples specifying the bounds of each dimension.
        D::Int64 # is the number of dimensions of the system, inferred from 'dimensions'.

    function System(particles,Hamiltonian,dimensions=fill((-10.,10.),3))

        for bounds in dimensions
            @assert (bounds[1] ≤ bounds[2]) "The given dimension bounds $(bounds[1]) to $(bounds[2]) are not valid; "*
                "input tuples in 'dimensions' such that the first number is lower than or equal to the second number."
        end

        for term in Hamiltonian , combination in term.groups , i in combination
            @assert (i ≤ length(particles)) "The given particle group index $i is not valid; "*
                "input only group indices in the energy terms of 'Hamiltonian' "*
                "smaller than or equal to the number of particle groups in 'particles'."
        end

        N = sum(group[1] for group in particles)
        D = length(dimensions)

        Hamiltonian = adjust_dimensions(Hamiltonian,D)

        new(particles,N,Hamiltonian,dimensions,D)
    end

end



HarmonicOscillator(particles,dimensions=fill((-10.,10.),3);strength=1.0,centre=[0.,0.,0.,]) =
    System(particles,HarmonicHamiltonian(strength,centre),dimensions)



; # suppresses inclusion output.
