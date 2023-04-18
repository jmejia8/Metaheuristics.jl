"""
    BitFlipMutation(;p = 1e-2)

Flip each bit with probability `p`.
"""
struct BitFlipMutation
    p::Float64
    rng
end

BitFlipMutation(;p = 1e-2, rng = default_rng_mh()) = BitFlipMutation(p, rng)

function mutation!(Q::AbstractMatrix{Bool}, parameters::BitFlipMutation)
    mask = rand(parameters.rng, size(Q)...) .< parameters.p
    Q[mask] = .!Q[mask]
    Q
end
