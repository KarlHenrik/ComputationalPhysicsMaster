# ENERGY 1-BODY TERMS
# to be used with the file ENERGY TERMS in the module QUANTUM DYNAMICS
# by Jonas Boym Flaten



abstract type Energy1Term <: EnergyTerm end # comprises various possible energy terms of the system Hamiltonian.



struct CustomEnergy1Term <: Energy1Term
    # represents arbitrary 1-body energy terms for the given particle groups, specified by involved properties and form.

    groups::Vector{Int64} # are the particle groups involved with the energy term.
    properties::Vector{ParticleProperty} # are the particle properties involved in the energy term.
    form::Function # is the mathematical form of the energy term, declared in the same order as the particle properties above.

    function CustomEnergy1Term(groups,properties,form)

        form_arguments = first(methods(form)).nargs-1 # is the number of arguments in the mathematical form of the energy term.
        @assert (form_arguments == length(properties)) "The declared number of arguments "*
            "does not match the declared number of properties; input as many elements in 'properties' as there are arguments in 'form'."

        return new(groups,properties,form)
    end
end

CustomEnergy1Term(groups;properties,form) = CustomEnergy1Term(groups,properties,form)
CustomEnergy1Term(group;properties,form) = CustomEnergy1Term([group],properties,form)
CustomEnergy1Term(groups...;properties,form) = CustomEnergy1Term([groups...],properties,form)

CustomEnergy1Terms(properties,form) = CustomEnergy1Term([0],properties,form)
    # returns a CustomEnergy1Term applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
CustomEnergy1Terms(;properties,form) = CustomEnergy1Terms(properties,form)



struct KineticTerm <: Energy1Term
    # represents one-body kinergy terms for the given particle groups, of the form p²/2m.

    groups::Vector{Int64} # are the particle groups involved with a kinetic term.

end

KineticTerm(group) = KineticTerm([group])
KineticTerm(groups...) = KineticTerm([groups...])

KineticTerms() = KineticTerm([0])
    # returns a KineticTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)



abstract type PotentialTerm <: Energy1Term end # comprises various possible one-body potergy terms of the system Hamiltonian.


struct HarmonicTerm <: PotentialTerm
    # represents one-body harmonic potergy terms for the given particle groups,
        # of the form m(ω²/2)(r-r₀)² with a strength ω²/2 and a center r₀.

    groups::Vector{Int64} # are the particle groups involved with the harmonic term.
    strength::Float64 # is the strength ω²/2 of the harmonic term.
    center::Vector{Float64} # is the center r₀ of the harmonic term.

    function HarmonicTerm(groups,strength=1.0,center=[0.,0.,0.])

        @assert (strength > 0) "The given harmonic strength $strength is invalid; input a positive number as 'strength'."

        return new(groups,strength,center)
    end
end

HarmonicTerm(groups::Vector{Int64};strength=1.0,center=[0.,0.,0.]) = HarmonicTerm(groups,strength,center)
HarmonicTerm(group;strength,center=[0.,0.,0.]) = HarmonicTerm([group],strength,center)
HarmonicTerm(groups...;strength,center=[0.,0.,0.]) = HarmonicTerm([groups...],strength,center)

HarmonicTerms(strength=1.0,center=[0.,0.,0.]) = HarmonicTerm([0],strength,center)
    # returns a HarmonicTerm applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
HarmonicTerms(;strength,center=[0.,0.,0.]) = HarmonicTerms(strength,center)

adjust_dimensions(term::HarmonicTerm,D::Int64) = HarmonicTerm(term.groups,term.strength,term.center[1:D])
    # adjusts the dimensions of the center of the HarmonicTerm to match the dimensions of the System.


struct Newton1Term <: PotentialTerm
    # represents one-body Newton potergy terms for the given particle groups,
        # of the form Am/√((r-r₀)²+a²) with a strength A, a center r₀ and a shielding a.

    groups::Vector{Int64} # are the particle groups involved with the Newton term.
    strength::Float64 # is the strength A of the Newton term.
    center::Vector{Float64} # is the center r₀ of the Newton term.
    shielding::Float64 # is the shielding a of the Newton term.

    function Newton1Term(groups,strength=1.0,center=[0.,0.,0.],shielding=0.01)

        @assert (strength > 0) "The given Newton strength $strength is invalid; input a positive number as 'strength'."

        return new(groups,strength,center,shielding)
    end
end

Newton1Term(groups::Vector{Int64};strength=1.0,center=[0.,0.,0.],shielding=0.01) = Newton1Term(groups,strength,center,shielding)
Newton1Term(group;strength,center=[0.,0.,0.],shielding=0.01) = Newton1Term([group],strength,center,shielding)
Newton1Term(groups...;strength,center=[0.,0.,0.],shielding=0.01) = Newton1Term([groups...],strength,center,shielding)

Newton1Terms(strength=1.0,center=[0.,0.,0.],shielding=0.01) = Newton1Term([0],strength,center,shielding)
    # returns a Newton1Term applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
Newton1Terms(;strength,center=[0.,0.,0.],shielding=0.01) = Newton1Terms(strength,center,shielding)

adjust_dimensions(term::Newton1Term,D::Int64) = Newton1Term(term.groups,term.strength,term.center[1:D],term.shielding)
    # adjusts the dimensions of the center of the Newton1Term to match the dimensions of the System.


struct Coulomb1Term <: PotentialTerm
    # represents one-body Coulomb potergy terms for the given particle groups,
        # of the form Aq/√((r-r₀)²+a²) with a strength A, a center r₀ and a shielding a.

    groups::Vector{Int64} # are the particle groups involved with the Coulomb term.
    strength::Float64 # is the strength A of the Coulomb term.
    center::Vector{Float64} # is the center r₀ of the Coulomb term.
    shielding::Float64 # is the shielding a of the Coulomb term.

    function Coulomb1Term(groups,strength=1.0,center=[0.,0.,0.],shielding=0.01)

        @assert (strength > 0) "The given Coulomb strength $strength is invalid; input a positive number as 'strength'."

        return new(groups,strength,center,shielding)
    end
end

Coulomb1Term(groups::Vector{Int64};strength=1.0,center=[0.,0.,0.],shielding=0.01) = Coulomb1Term(groups,strength,center,shielding)
Coulomb1Term(group;strength,center=[0.,0.,0.],shielding=0.01) = Coulomb1Term([group],strength,center,shielding)
Coulomb1Term(groups...;strength,center=[0.,0.,0.],shielding=0.01) = Coulomb1Term([groups...],strength,center,shielding)

Coulomb1Terms(strength=1.0,center=[0.,0.,0.],shielding=0.01) = Coulomb1Term([0],strength,center,shielding)
    # returns a Coulomb1Term applying to all particles in the system. (groups = [0] encodes that the term applies to all particles.)
Coulomb1Terms(;strength,center=[0.,0.,0.],shielding=0.01) = Coulomb1Terms(strength,center,shielding)

adjust_dimensions(term::Coulomb1Term,D::Int64) = Coloumb1Term(term.groups,term.strength,term.center[1:D],term.shielding)
    # adjusts the dimensions of the center of the Newton1Term to match the dimensions of the System.



; # suppresses inclusion output.
