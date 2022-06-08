abstract type Interaction end

struct CalogeroSutherland <: Interaction
    beta_kinda::Float64
    function CalogeroSutherland(beta)
        return new(beta * (beta - 1))
    end
end

function interaction_over_grid!(interaction, x1, grid, V::CalogeroSutherland)
    interaction .= V.beta_kinda ./ (grid .- x1).^2
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