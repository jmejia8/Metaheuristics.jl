"""
    SBX(;η, p, bounds)

Simulated Binomial Crossover.
"""
mutable struct SBX
    η::Float64
    p::Float64
    bounds
    rng
end

SBX(;η=15, p=0.9, bounds=zeros(0,0), rng = default_rng_mh()) = SBX(η, p, bounds, rng)


function crossover(population, parameters::SBX)
    isempty(population) && return zeros(0,0)
    Q = positions(population)
    bounds = parameters.bounds
    rng = parameters.rng
    for i in 1:2:length(population)-1
        p1 = get_position(population[i])
        p2 = get_position(population[i+1])
        c1, c2 = SBX_crossover(p1, p2, bounds, parameters.η, parameters.p, rng)
        Q[i,:] = c1
        Q[i+1,:] = c2
    end
    Q
end

function gen_β(β, η, D, R)
    α = 2 .- β .^ (-  η - 1 )
    mask = R .<= 1 ./ α
    s = 1 / (η + 1)
    βq = [ mask[i] ?  (R[i] * α[i])^s : (1 / (2- R[i]*α[i]))^s for i in 1:D]
    βq
end

"""
    SBX_crossover(vector1, vector2, bounds, η=15, p_variable = 0.9)

Simulated binomial crossover for given two `Vectors{Real}`.
"""
function SBX_crossover(
        vector1,
        vector2,
        bounds::BoxConstrainedSpace,
        η=15,
        p_variable = 0.9,
        rng = default_rng_mh()
    )

    xu = bounds.ub
    xl = bounds.lb
    D = getdim(bounds)

    do_crossover = ones(Bool, D)
    do_crossover[rand(rng, D) .> p_variable] .= false
    do_crossover[ abs.( vector2 - vector1 ) .<= eps() ] .= false

    y1 = min.( vector1, vector2 )
    y2 = max.( vector1, vector2 )
    Δ = max.(eps(), y2 - y1)


    R = rand(rng, D)

    β = @. 1 + (2* (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c1 = @. 0.5*(y1 + y2 -  βq*Δ)

    β = @. 1 + (2* (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c2 = @. 0.5*(y1 + y2 +  βq*Δ)

    # swap
    mask = rand(rng, Bool, D)
    cc = copy(c1)
    c1[mask] = c2[mask]
    c2[mask] = cc[mask]

    cc1 = copy(vector1)
    cc1[do_crossover] = _to_int_if_necessary(eltype(cc1), c1[do_crossover] )
    cc2 = copy(vector2)
    cc2[do_crossover] = _to_int_if_necessary(eltype(cc2), c2[do_crossover] )


    reset_to_violated_bounds!(cc1, bounds)
    reset_to_violated_bounds!(cc2, bounds)

    return cc1, cc2
end

function SBX_crossover(
        vector1,
        vector2,
        bounds::AbstractMatrix,
        η=15,
        p_variable = 0.9,
        rng = default_rng_mh()
    )
    b = BoxConstrainedSpace(lb = bounds[1,:], ub = bounds[2,:])
    SBX_crossover(vector1, vector2, b, η, p_variable, rng)
end

