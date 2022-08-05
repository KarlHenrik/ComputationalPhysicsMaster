abstract type Interaction end

struct CalogeroSutherland <: Interaction
    ββ_1::Float64
    function CalogeroSutherland(beta)
        ββ_1 = beta * (beta - 1)
        return new(ββ_1)
    end
end

function interaction_over_grid!(interaction, x1, grid, V::CalogeroSutherland)
    interaction .= V.ββ_1 ./ ((grid .- x1).^2 .+ 0.1)
end

struct ShieldedCoulomb <: Interaction
    shielding::Float64
end

function interaction_over_grid!(interaction, x1, grid, V::ShieldedCoulomb)
    interaction .= 1 ./ sqrt.( (grid .- x1).^2 .+ V.shielding.^2 )
    return interaction
end

struct Coulomb <: Interaction end

function interaction_over_grid!(interaction, x1, grid, V::Coulomb)
    interaction .= 1 ./ sqrt.( (grid .- x1).^2 )
    return interaction
end

struct NonInteracting <: Interaction end