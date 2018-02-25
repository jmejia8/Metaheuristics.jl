# struct Individual
#     x::Vector
#     f::Real
#     g::Vector
#     h::Vector
#     νVal::Real
# end

# function reducePop(Population, N, N_new)
#     fitness = zeros(N)
#     for i = 1:N
#         fitness[i] = Population[i].f
#     end

#     fitness_order = sort(fitness)

#     P = Array{Individual, 1}([])

#     for i = 1:N_new
#         item = findfirst(x -> x == fitness_order[i], fitness)
#         push!(P, Population[item])
#     end

#     return P

# end

# function getWorst(Population::Array{Individual, 1})
#     ν_max = Population[1].νVal
#     f_max = Population[1].f
#     j = 1

#     for i = 2:length(Population)
#         if ν_max > Population[i].νVal 
#             ν_max = Population[i].νVal
#             f_max = Population[i].f
#             j = i
#         elseif ν_max == Population[i].νVal && f_max < Population[i].f
#             ν_max = Population[i].νVal
#             f_max = Population[i].f
#             j = i
#         end
#     end

#     return Population[j]
# end

###################################################
#      Solutions and population functions
#          for Matrix representation
###################################################
function correctSol(y::Vector, a::Vector, b::Vector)
    # Correct solution

    for i = 1:length(y)
        if !( a[i] <= y[i] <= b[i] )
            y[i] = a[i] + (b[i] - a[i])*rand()
        end
    end
    
    return y
end

function correctPop(P, a, b)
    # Correct population
    # a, b should be D × 1
    N, D = size(P)

    for i = 1:N
        for j = 1:D
            if !(a[j] <= P[i, j] <= b[j])
                P[i, j] = a[j] - (b[j] - a[j])*rand()
            end
        end
    end

    return P
end

function initializePop(N, D, a, b)
    # a, b should be D × 1
    return a' .* ones(N, D) + (b - a)' .* rand(N, D)
end

function initializeSol(D, a, b)
    # a, b should be D × 1
    return a + (b - a) .* rand(D)
end

###################################################
#             Fitness functions
#          for vector representation
###################################################
function evaluatePop(X::Matrix, fobj::Function, N::Int)

    f = zeros(N)
    for i = 1:N
        f[i] = fobj(X[i,:])
    end

    return f
end

function getBest(fitness::Vector, searchType::Symbol = :minimize)
    if searchType == :minimize
        best_X = indmin(fitness) # minimization.
        best = fitness[best_X] 
    else
        best_X = indmax(fitness) # maximization.
        best = fitness[best_X] 
    end

    return best_X, best
end