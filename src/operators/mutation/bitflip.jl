"""
    BitFlipMutation(;p = 1e-2)

Flip each bit with probability `p`.
"""
struct BitFlipMutation{R}
    p::Float64
    rng::R
end

BitFlipMutation(p, rng) = BitFlipMutation{typeof(rng)}(p, rng)

BitFlipMutation(; p = 1.0e-2, rng = default_rng_mh()) = BitFlipMutation(p, rng)

function mutation!(Q::AbstractMatrix{Bool}, parameters::BitFlipMutation)
    # one bulk draw keeps seeded results unchanged
    R = rand(parameters.rng, size(Q)...)
    p = parameters.p
    @inbounds for k in eachindex(Q, R)
        Q[k] = xor(Q[k], R[k] < p) # branchless flip
    end
    Q
end
