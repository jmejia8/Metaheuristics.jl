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


function getU(P::Array{xfgh_indiv}, K::Int, I::Vector{Int}, i::Int, N::Int, feasible_solutions)
    # at least half of the population is feasible to generate random centers
    if length(feasible_solutions) >= 0.5N || length(feasible_solutions) == 0
        return getU(P, K, I, i, N)
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

    return P[U_ids]
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

