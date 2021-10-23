"""
    GA_reproduction(pa::AbstractVector{T},
                    pb::AbstractVector{T},
                    bounds;
                    η_cr = 20,
                    η_m  = 15,
                    p_cr = 0.9,
                    p_m  = 0.1)

Crate two solutions by applying SBX to parents pa and pb and polynomial mutation to
offspring. Return two vectors.
"""
function GA_reproduction(pa::AbstractVector{T},
        pb::AbstractVector{T},
        bounds;
        η_cr = 20,
        η_m  = 15,
        p_cr = 0.9,
        p_m  = 0.1
    ) where T <: AbstractFloat


    # crossover
    c1, c2 = SBX_crossover(pa, pb, bounds, η_cr, p_cr)

    # mutation
    polynomial_mutation!(c1, bounds,η_m, p_m)
    polynomial_mutation!(c2, bounds,η_m, p_m)

    # rapair solutions if necesary
    reset_to_violated_bounds!(c1, bounds)
    reset_to_violated_bounds!(c2, bounds)

    return c1, c2
end


"""
    GA_reproduction_half(pa::AbstractVector{T},
                    pb::AbstractVector{T},
                    bounds;
                    η_cr = 20,
                    η_m  = 15,
                    p_cr = 0.9,
                    p_m  = 0.1)

Same that `GA_reproduction` but only returns one offspring.
"""
function GA_reproduction_half(pa::AbstractVector{T},
        pb::AbstractVector{T},
        bounds;
        η_cr = 20,
        η_m  = 15,
        p_cr = 0.9,
        p_m  = 0.1
    ) where T <: AbstractFloat


    # crossover
    _, c2 = SBX_crossover(pa, pb, bounds, η_cr, p_cr)

    # mutation
    polynomial_mutation!(c, bounds,η_m, p_m)

    # rapair solution if necesary
    reset_to_violated_bounds!(c, bounds)

    return c
end
