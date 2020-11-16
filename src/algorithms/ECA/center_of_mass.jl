function center(U, mass)
    d = length(U[1].x)

    c = zeros(Float64, d)

    for i = 1:length(mass)
        c += mass[i] .* U[i].x
    end

    return c / sum(mass)
end

function center(U::Array, searchType::Symbol; ε = 0.0)
    n = length(U)

    mass = getMass(U, searchType; ε = ε)

    return center(U, mass),
    getWorstInd(U, searchType, is_better_eca),
    getBestInd(U, searchType, is_better_eca)
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

