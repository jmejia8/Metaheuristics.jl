################################################
# Unconstrained
################################################
function fitnessToMass(fitness::Vector{Float64})
    m = minimum(fitness)

    if m < 0
        fitness = 2 * abs(m) .+ fitness
    end
    
    fitness = 2 * maximum(fitness) .- fitness

    return fitness
end

"""
    getMass(U, max_fitness, epsilon)
    get the mass in a vector `m[i] = f[i] + 2*max_fitness*sum_vio[i] `, tolerance
    in equality constraints are given by `epsilon`
"""
function getMass(U::Array{xf_indiv,1})
    n, d = length(U), length(U[1].x)

    fitness = zeros(Float64, n)

    for i = 1:n
        fitness[i] = U[i].f
    end

    return fitnessToMass(fitness)
end


function center(U, mass)
    d = length(U[1].x)

    c = zeros(Float64, d)

    for i = 1:length(mass)
        c += mass[i] .* U[i].x
    end

    return c / sum(mass)
end

function center(U::Array)
    n = length(U)

    mass = getMass(U)

    return center(U, mass),
    argworst(U), # worst
    argbest(U)  # best
end

function getU(P::Array, K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N - K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j.+1]
    end

    return P[U_ids]
end


function getU_ids(K::Int, I::Vector{Int}, i::Int, N::Int)
    if i <= N - K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j.+1]
    end

    return U_ids
end



function getU_ids(K::Int, I::Vector{Int}, i::Int, N::Int, feasible_solutions)
    # at least half of the population is feasible to generate random centers
    if length(feasible_solutions) >= 0.5N || length(feasible_solutions) == 0
        return getU_ids(K, I, i, N)
    end
    
    # center of mass is generated with at least one feasible solution
    K -= 1
    if i <= N - K
        U_ids = I[i:K+i]
    else
        j = (i:K+i) .% N
        U_ids = I[j.+1]
    end
    
    push!(U_ids, rand(feasible_solutions))

    return U_ids
end

function crossover(
    x::Vector{Float64},
    y::Vector{Float64},
    p_cr::Vector{Float64},
)
    D = length(x)
    tmp2 = zeros(D)
    for j = 1:D
        if rand() < p_cr[j]
            y[j] = x[j]
            tmp2[j] += 1
        end
    end

    return y, tmp2
end


"""
    ECA_operator(population, K, η_max)

Compute a solution using ECA variation operator, `K` is the number of solutions used to
calculate the center of mass and `η_max` is the maximum stepsize.
"""
function ECA_operator(
        population::AbstractArray{xf_indiv}, K, η_max;
        i = rand(1:length(population)),
        U = rand(population, K),
        bounds = nothing
    )

    x = get_position(population[i])

    # generate center of mass
    c, u_worst, u_best = center(U)

    # stepsize
    η = η_max * rand()

    # u: worst element in U
    u = get_position(U[u_worst])

    y = x .+ η .* (c .- u)
    if isnothing(bounds)
        return y
    end

    evo_boundary_repairer!(y, c, bounds)

    return y
end

