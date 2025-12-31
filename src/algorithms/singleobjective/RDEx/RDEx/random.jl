"""
Random number generation utilities for RDEx optimizer.

This module encapsulates all random number generation to avoid global state.
Uses a single RNG for all random operations, which simplifies the code and
reduces memory usage while maintaining reproducibility.
"""
module RandomUtils

using Random
using Distributions

"""
    RNGState

Encapsulates a single random number generator used by the optimizer.
This avoids global state and allows for proper seeding and reproducibility.
"""
struct RNGState
    rng::MersenneTwister
end

"""
    RNGState(seed::Integer)

Create a new RNGState with the given seed.
"""
function RNGState(seed::Integer)
    seed_u = UInt32(seed)
    return RNGState(MersenneTwister(seed_u))
end

"""
    RNGState()

Create a new RNGState with a random seed.
"""
function RNGState()
    return RNGState(MersenneTwister())
end

"""
    int_random(rng::RNGState, target::Int) -> Int

Generate a random integer in [0, target-1] if target > 0, else 0.
"""
function int_random(rng::RNGState, target::Int)::Int
    target <= 0 && return 0
    return rand(rng.rng, 0:(target-1))
end

"""
    random_val(rng::RNGState, minimal::Float64, maximal::Float64) -> Float64

Generate a random value uniformly distributed in [minimal, maximal).
"""
function random_val(rng::RNGState, minimal::Float64, maximal::Float64)::Float64
    return rand(rng.rng) * (maximal - minimal) + minimal
end

"""
    norm_rand(rng::RNGState, mu::Float64, sigma::Float64) -> Float64

Generate a random value from a normal distribution N(mu, sigma).
"""
function norm_rand(rng::RNGState, mu::Float64, sigma::Float64)::Float64
    return randn(rng.rng) * sigma + mu
end

"""
    cachy_rand(rng::RNGState, mu::Float64, sigma::Float64) -> Float64

Generate a random value from a Cauchy distribution with location mu and scale sigma.
"""
function cachy_rand(rng::RNGState, mu::Float64, sigma::Float64)::Float64
    return rand(rng.rng, Cauchy(mu, sigma))
end

end # module RandomUtils

