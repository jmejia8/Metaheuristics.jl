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


end # module RandomUtils

