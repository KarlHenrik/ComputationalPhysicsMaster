# ENERGY 2-BODY TERMS
# to be used with the file ENERGY TERMS in the module QUANTUM DYNAMICS
# by Jonas Boym Flaten



abstract type Energy2Term <: EnergyTerm end # comprises various possible energy terms of the system Hamiltonian.



struct CustomEnergy2Term <: Energy2Term
    # represents arbitrary 2-body energy terms for the given particle groups, specified by involved properties and form.

    groups::Vector{Vector{Int64}} # are the combinations of particle groups involved with the energy term.
    properties::Vector{Vector{ParticleProperty}} # are the particle properties involved in the energy term, sorted by body.
    form::Function # is the mathematical form of the energy term, declared in the same order as the particle properties above.

    function CustomEnergy2Term(groups,properties,form)

        form_arguments = first(methods(form)).nargs-1 # is the number of arguments in the mathematical form of the energy term.
        @assert (form_arguments == sum([length(body) for body in properties])) "The declared number of arguments "*
            "does not match the declared number of properties; input as many elements in 'properties' as there are arguments in 'form'."

        return new(groups,properties,form)
    end
end

CustomEnergy2Term(groups;properties,form) = CustomEnergy2Term(groups,properties,form)

function CustomEnergy2Term(within_groups::Vector{Int64},between_groups::Vector{Int64},properties,form)
    # returns a CustomEnergy2Term which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.

    bodies = length(properties) # is the number of bodies involved in the energy term.
    groups = [within(within_groups)...,between(between_groups;bodies=bodies)...]

    return CustomEnergy2Term(groups,properties,form)
end

CustomEnergy2Term(within_groups,between_groups;properties,form) = CustomEnergy2Term(within_groups,between_groups,properties,form)

CustomEnergy2Terms(properties,form) = CustomEnergy2Term([[0,0]],properties,form)
    # returns a CustomEnergy2Term applying to all particles in the system. (groups = [[0,0]] encodes that the term applies to all particles.)
CustomEnergy2Terms(;properties,form) = CustomEnergy2Terms(properties,form)



struct Newton2Term <: Energy2Term
    # represents two-body Newton interaction terms for the given particle groups,
        # of the form Am₁m₂/√((r₁-r₂)²+a²) with a strength A and a shielding a.

    groups::Vector{Vector{Int64}} # are the combinations of particle groups involved with Newton interaction.
    strength::Float64 # is the strength A of the Newton term.
    shielding::Float64 # is the shielding a of the Newton term.

    function Newton2Term(groups,strength,shielding)

        @assert (strength > 0) "The given Newton interaction strength $strength is invalid; input a positive number as 'strength'."

        return new(groups,strength,shielding)
    end
end

Newton2Term(groups::Vector{Tuple{Int64,Int64}};strength=1.0,shielding=0.01) = Newton2Term(groups,strength,shielding)

function Newton2Term(within_groups::Vector{Int64},between_groups::Vector{Int64},strength::Float64=1.0,shielding::Float64=0.01)
    # returns a Newton2Term which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.

    groups = [within(within_groups)...,between(between_groups;bodies=2)...]

    return Newton2Term(groups,strength,shielding)
end

Newton2Term(within_groups,between_groups;strength,shielding=0.01) = Newton2Term(within_groups,between_groups,strength,shielding)

Newton2Terms(strength=1.0,shielding=0.01) = Newton2Term([[0,0]],strength,shielding)
    # returns a Newton2Term applying to all particles in the system. (groups = [[0,0]] encodes that the term applies to all particles.)
Newton2Terms(;strength,shielding=0.01) = Newton2Terms(strength,shielding)



struct Coulomb2Term <: Energy2Term
    # represents two-body Coulomb interaction terms for the given particle groups,
        # of the form Aq₁q₂/√((r₁-r₂)²+a²) with a strength A and a shielding a.

    groups::Vector{Vector{Int64}} # are the combinations of particle groups involved with Coulomb interaction.
    strength::Float64 # is the strength A of the Coulomb term.
    shielding::Float64 # is the shielding a of the Coulomb term.

    function Coulomb2Term(groups,strength=1.0,shielding=0.01)

        @assert (strength > 0) "The given Coulomb interaction strength $strength is invalid; input a positive number as 'strength'."

        return new(groups,strength,shielding)
    end
end

Coulomb2Term(groups::Vector{Tuple{Int64,Int64}};strength,shielding=0.01) = Coulomb2Term(groups,strength,shielding)

function Coulomb2Term(within_groups::Vector{Int64},between_groups::Vector{Int64},strength::Float64=1.0,shielding::Float64=0.01)
    # returns a Coulomb2Term which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.

    groups = [within(within_groups)...,between(between_groups;bodies=2)...]

    return Coulomb2Term(groups,strength,shielding)
end

Coulomb2Term(within_groups,between_groups;strength,shielding=0.01) = Coulomb2Term(within_groups,between_groups,strength,shielding)

Coulomb2Terms(strength=1.0,shielding=0.01) = Coulomb2Term([[0,0]],strength,shielding)
    # returns a Coulomb2Term applying to all particles in the system. (groups = [[0,0]] encodes that the term applies to all particles.)
Coulomb2Terms(;strength,shielding=0.01) = Coulomb2Terms(strength,shielding)



; # suppresses inclusion output.
