"""
Mutation strategies for RDEx optimizer.

Implements both standard and elite-based mutation strategies.
"""
module MutationUtils

using Random
using Distributions
import ..Constants

"""
    EliteBasedOrdering

Stores ordered positions for elite-based mutation.
"""
struct EliteBasedOrdering
    best::Vector{Float64}
    medium::Vector{Float64}
    worst::Vector{Float64}
end

"""
    elite_based_order(
        pos1::Vector{Float64}, pos1_fit::Float64,
        pos3::Vector{Float64}, pos3_fit::Float64,
        pos4::Vector{Float64}, pos4_fit::Float64
    ) -> EliteBasedOrdering

Order three positions by fitness for elite-based mutation.
Returns the archive ordering.
"""
function elite_based_order(
    pos1::Vector{Float64}, pos1_fit::Float64,
    pos3::Vector{Float64}, pos3_fit::Float64,
    pos4::Vector{Float64}, pos4_fit::Float64
)::EliteBasedOrdering
    return _order_three(pos1, pos1_fit, pos3, pos3_fit, pos4, pos4_fit)
end

function _order_three(
    pos1::Vector{Float64}, pos1_fit::Float64,
    pos3::Vector{Float64}, pos3_fit::Float64,
    pos4::Vector{Float64}, pos4_fit::Float64
)::EliteBasedOrdering
    if pos1_fit <= pos3_fit && pos1_fit <= pos4_fit
        best = pos1
        if pos3_fit <= pos4_fit
            return EliteBasedOrdering(best, pos3, pos4)
        else
            return EliteBasedOrdering(best, pos4, pos3)
        end
    elseif pos3_fit <= pos1_fit && pos3_fit <= pos4_fit
        best = pos3
        if pos1_fit <= pos4_fit
            return EliteBasedOrdering(best, pos1, pos4)
        else
            return EliteBasedOrdering(best, pos4, pos1)
        end
    else  # pos4_fit <= pos1_fit && pos4_fit <= pos3_fit
        best = pos4
        if pos1_fit <= pos3_fit
            return EliteBasedOrdering(best, pos1, pos3)
        else
            return EliteBasedOrdering(best, pos3, pos1)
        end
    end
end

"""
    create_elite_based_trial(
        x_chosen::Vector{Float64},
        ordering::EliteBasedOrdering,
        F::Float64,
        Cr::Float64,
        lower_bounds::Vector{Float64},
        upper_bounds::Vector{Float64},
        rng::AbstractRNG
    ) -> Tuple{Vector{Float64}, Float64}

Create trial vector using elite-based mutation strategy.
Returns (trial_vector, actual_crossover_rate).
"""
function create_elite_based_trial(
    x_chosen::Vector{Float64},
    ordering::EliteBasedOrdering,
    F::Float64,
    Cr::Float64,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    rng::AbstractRNG
)::Tuple{Vector{Float64}, Float64}
    D = length(x_chosen)
    trial = zeros(D)
    will_crossover = rand(rng, 0:(D-1))
    perturbation = rand(rng) < 0.4
    actual_cr = 0.0
    
    for j = 1:D
        if rand(rng) < Cr || will_crossover == j - 1
            trial[j] = x_chosen[j] +
                      F * (ordering.best[j] - x_chosen[j]) +
                      F * (ordering.medium[j] - ordering.worst[j])
            
            # Boundary handling
            if trial[j] < lower_bounds[j]
                trial[j] = rand(rng) * (upper_bounds[j] - lower_bounds[j]) + lower_bounds[j]
            elseif trial[j] > upper_bounds[j]
                trial[j] = rand(rng) * (upper_bounds[j] - lower_bounds[j]) + lower_bounds[j]
            end
            actual_cr += 1.0
        else
            trial[j] = perturbation ? 
                rand(rng, Cauchy(x_chosen[j], 0.1)) :
                x_chosen[j]
        end
    end
    
    return (trial, actual_cr / D)
end

"""
    create_standard_trial(
        x_chosen::Vector{Float64},
        x_prand::Vector{Float64},
        x_rand1::Vector{Float64},
        x_rand2::Vector{Float64},
        F::Float64,
        Cr::Float64,
        lower_bounds::Vector{Float64},
        upper_bounds::Vector{Float64},
        rng::AbstractRNG
    ) -> Tuple{Vector{Float64}, Float64}

Create trial vector using standard mutation strategy.
Returns (trial_vector, actual_crossover_rate).
"""
function create_standard_trial(
    x_chosen::Vector{Float64},
    x_prand::Vector{Float64},
    x_rand1::Vector{Float64},
    x_rand2::Vector{Float64},
    F::Float64,
    Cr::Float64,
    lower_bounds::Vector{Float64},
    upper_bounds::Vector{Float64},
    rng::AbstractRNG
)::Tuple{Vector{Float64}, Float64}
    D = length(x_chosen)
    trial = zeros(D)
    will_crossover = rand(rng, 0:(D-1))
    perturbation = rand(rng) < 0.4
    actual_cr = 0.0
    
    for j = 1:D
        if rand(rng) < Cr || will_crossover == j - 1
            trial[j] = x_chosen[j] +
                      F * (x_prand[j] - x_chosen[j]) +
                      F * (x_rand1[j] - x_rand2[j])
            
            # Boundary handling
            if trial[j] < lower_bounds[j]
                trial[j] = rand(rng) * (upper_bounds[j] - lower_bounds[j]) + lower_bounds[j]
            elseif trial[j] > upper_bounds[j]
                trial[j] = rand(rng) * (upper_bounds[j] - lower_bounds[j]) + lower_bounds[j]
            end
            actual_cr += 1.0
        else
            trial[j] = perturbation ?
                rand(rng, Cauchy(x_chosen[j], 0.1)) :
                x_chosen[j]
        end
    end
    
    return (trial, actual_cr / D)
end

end # module MutationUtils

