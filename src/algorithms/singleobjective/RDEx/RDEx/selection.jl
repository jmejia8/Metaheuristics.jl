"""
Selection and sampling utilities for RDEx optimizer.
"""
module SelectionUtils

using Random

"""
    weighted_sample(rng::AbstractRNG, weights::Vector{Float64}) -> Int

Perform weighted random sampling from a discrete distribution.

# Arguments
- `rng`: Random number generator
- `weights`: Vector of weights (non-negative, can sum to any value)

# Returns
0-based index of the selected element.
"""
function weighted_sample(rng::AbstractRNG, weights::Vector{Float64})::Int
    total = sum(weights)
    if total == 0.0 || isempty(weights)
        return rand(rng, 0:(length(weights)-1))
    end
    
    r = rand(rng) * total
    cumsum = 0.0
    for (i, w) in enumerate(weights)
        cumsum += w
        if r <= cumsum
            return i - 1  # 0-indexed
        end
    end
    return length(weights) - 1
end

end # module SelectionUtils

