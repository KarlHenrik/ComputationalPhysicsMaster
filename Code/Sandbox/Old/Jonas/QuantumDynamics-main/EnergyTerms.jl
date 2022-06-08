# ENERGY TERMS
# to be used with the module QUANTUM DYNAMICS
# by Jonas Boym Flaten



function within(groups::Vector{Int64};bodies::Int64=2)
    # returns a vector of tuples of the given particle group indices with the given number of bodies.

    return [fill(group,bodies) for group in groups]
end

within(groups...;bodies=2) = within([groups...];bodies=bodies)


function between(groups::Vector{Int64};bodies::Int64=2)
    # returns a vector of all combinations from the given particle groups with the given number of bodies.

    @assert (bodies ≤ length(groups)) "The number of involved bodies is larger than the number of given particle groups; "*
        "input at least as many indices in 'groups' as 'bodies'."

    if (bodies == 1)
        return [[group] for group in groups]
    else

        combinations = Vector{Int64}[]
        for i in 1:(length(groups)-(bodies-1))
            for combination in between(groups[(i+1):end];bodies=(bodies-1))
                push!(combinations,[groups[i],combination...])
            end
        end

        return combinations
    end
end

between(groups...;bodies=2) = between([groups...];bodies=bodies)
bethreen(groups) = between(groups;bodies=3)
bethreen(groups...) = bethreen([groups...])


abstract type EnergyTerm end # comprises various possible energy terms of the system Hamiltonian.

adjust_dimensions(term::EnergyTerm,D::Int64) = term

function adjust_dimensions(Hamiltonian::Vector{Term},D::Int64) where Term <: EnergyTerm

    return adjust_dimensions.(Hamiltonian,D)
end


@enum ParticleProperty m q s r p l # encodes the particle properties involved in energy terms.


include("Energy1Terms.jl") # includes the various possible 1-body energy terms of the system Hamiltonian defined in EnergyTerms.jl.
include("Energy2Terms.jl") # includes the various possible 2-body energy terms of the system Hamiltonian defined in EnergyTerms.jl.



struct CustomEnergyTerm <: EnergyTerm
    # represents arbitrary energy terms for the given particle groups, specified by involved properties and form.

    groups::Vector{Vector{Int64}} # are the combinations of particle groups involved with the energy term.
    properties::Vector{Vector{ParticleProperty}} # are the particle properties involved in the energy term, sorted by body.
    form::Function # is the mathematical form of the energy term, declared in the same order as the particle properties above.

    function CustomEnergyTerm(groups,properties,form)

        bodies = length(properties) # is the number of bodies involved in the energy term.
        for combination in groups
            @assert (length(combination) == bodies) "The number of involved bodies does not match the given particle group combinations; "*
                "input tuples in 'groups' whose sizes match the number of tuples in 'properties'."
        end

        form_arguments = first(methods(form)).nargs-1 # is the number of arguments in the mathematical form of the energy term.
        @assert (form_arguments == sum([length(body) for body in properties])) "The declared number of arguments "*
            "does not match the declared number of properties; input as many elements in 'properties' as there are arguments in 'form'."

        return new(groups,properties,form)
    end
end

CustomEnergyTerm(groups;properties,form) = CustomEnergyTerm(groups,properties,form)

function CustomEnergyTerm(within_groups::Vector{Int64},between_groups::Vector{Int64},properties,form)
    # returns a CustomEnergyTerm which applies within the particle groups given in 'within_groups'
        # and between the particle groups given in 'between_groups'.

    bodies = length(properties) # is the number of bodies involved in the energy term.
    groups = [within(within_groups)...,between(between_groups;bodies=bodies)...]

    return CustomEnergyTerm(groups,properties,form)
end

CustomEnergyTerm(within_groups,between_groups;properties,form) = CustomEnergyTerm(within_groups,between_groups,properties,form)

CustomEnergyTerms(properties,form) = CustomEnergyTerm([fill(0,length(properties))],properties,form)
    # returns a CustomEnergyTerm applying to all particles in the system. (groups = [[0,...,0]] encodes that the term applies to all particles.)
CustomEnergyTerms(;properties,form) = CustomEnergyTerms(properties,form)



FreeHamiltonian() = [KineticTerms()]

HarmonicHamiltonian(strength=1.0,center=[0.,0.,0.]) = [KineticTerms(),HarmonicTerms(strength,center)]
HarmonicHamiltonian(;strength=1.0,center=[0.,0.,0.]) = HarmonicHamiltonian(strength,center)

NewtonHamiltonian(strength=1.0,shielding=0.01) = [KineticTerms(),Newton2Terms(strength,shielding)]
NewtonHamiltonian(;strength=1.0,shielding=0.01) = NewtonHamiltonian(strength,shielding)

CoulombHamiltonian(strength=1.0,shielding=0.01) = [KineticTerms(),Coulomb2Terms(strength,shielding)]
CoulombHamiltonian(strength=1.0,shielding=0.01) = CoulombHamiltonian(strength,shielding)



; # suppresses inclusion output.
