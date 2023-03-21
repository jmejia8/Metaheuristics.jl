"""
    UniformCrossover(;p = 0.5)

Uniform crossover aka Binomial crossover suitable for binary representation.
"""
struct UniformCrossover
    p::Float64
    UniformCrossover(;p=0.5) = new(p)
end

function crossover(population, parameters::UniformCrossover)
    n = length(population) ÷ 2
    offspring_A = positions(population[1:n])
    offspring_B = positions(population[n+1:2n])
    mask = rand(size(offspring_A)...) .<= parameters.p
    tmp = copy(offspring_A[mask])
    offspring_A[mask] = offspring_B[mask]
    offspring_B[mask] = tmp
    [offspring_A; offspring_B]
end

"""
    OrderCrossover()
Order crossover for representation where order is important. Suitable for permutation
representation.
"""
struct OrderCrossover end

function crossover(population, parameters::OrderCrossover)
    O = positions(population)
    N, D = size(O)
    s = rand(1:D, N) # slash points
    for i = 1:2:N
        PA = O[i, :]   # parent A
        PB = O[i+1, :] # parent B
        O[i,  s[i]+1:D] = setdiff(PB, PA[1:s[i]]);
        O[i+1,s[i]+1:D] = setdiff(PA, PB[1:s[i]]);
    end
    O
end

##########################################################################

function gen_β(β, η, D, R)
    α = 2.0 .- β .^ (-  η - 1.0 )
    mask = R .<= 1.0 ./ α
    s = 1.0 / (η + 1.0)
    βq = [ mask[i] ?  (R[i] * α[i])^s : (1.0 / (2.0 - R[i]*α[i]))^s for i in 1:D]
    βq
end

"""
    SBX_crossover(vector1, vector2, bounds, η=15, p_variable = 0.9)

Simulated binomial crossover for given two `Vectors{Real}`.
"""
function SBX_crossover(vector1, vector2, bounds::Bounds, η=15, p_variable = 0.9)
    xu = bounds.ub
    xl = bounds.lb
    D = getdim(bounds)

    do_crossover = ones(Bool, D)
    do_crossover[rand(D) .> p_variable] .= false
    do_crossover[ abs.( vector2 - vector1 ) .<= eps() ] .= false

    y1 = min.( vector1, vector2 )
    y2 = max.( vector1, vector2 )
    Δ = max.(eps(), y2 - y1)


    R = rand(D)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c1 = @. 0.5*(y1 + y2 -  βq*Δ)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c2 = @. 0.5*(y1 + y2 +  βq*Δ)

    # swap
    mask = rand(Bool, D)
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

"""
    SBX(;η, p, bounds)

Simulated Binomial Crossover.
"""
mutable struct SBX
    η::Float64
    p::Float64
    bounds
    SBX(;η = 15, p = 0.9, bounds = zeros(0,0)) = new(η, p, bounds)
end

function crossover(population, parameters::SBX)
    isempty(population) && return zeros(0,0)
    Q = positions(population)
    bounds = parameters.bounds
    for i in 1:2:length(population)-1
        p1 = get_position(population[i])
        p2 = get_position(population[i+1])
        c1, c2 = SBX_crossover(p1, p2, bounds, parameters.η, parameters.p)
        Q[i,:] = c1
        Q[i+1,:] = c2
    end
    Q
end



"""
    DE_crossover(x, u, CR)

Binomial crossover between x and u for Differential Evolution with probability CR, i.e.,
`v[j] = u[j]` if `rand() < CR`, otherwise `v[j] = x[j]`. Return `v`.
"""
function DE_crossover(x, u, CR)
    D = length(x)
    # binomial crossover
    v = zeros(D)
    j_rand = rand(1:D)

    # binomial crossover
    for j = 1:D
        if rand() < CR || j == j_rand
            v[j] = u[j]
        else
            v[j] = x[j]
        end
    end

    return v
end

