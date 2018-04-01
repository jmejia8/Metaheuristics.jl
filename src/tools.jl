struct Individual_
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
end

###################################################
#      Solutions and population functions
#          for Matrix representation
###################################################
function correctSol(y::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64})
    # Correct solution

    for i = 1:length(y)
        if !( a[i] <= y[i] <= b[i] )
            y[i] = a[i] + (b[i] - a[i])*rand()
        end
    end
    
    return y
end


function correct(y::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64}, crrect::Bool=false)
    if crrect
        return correctSol(y, a, b)
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


function initializePop(N::Int, D::Int, a::Vector{Float64}, b::Vector{Float64}, initType::Symbol=:uniform)
    # a, b should be D × 1

    if initType == :cheb
        chebPts(x, a, b) = 0.5*(a + b) + 0.5*(b-a)*cos.( x )
        X = zeros(N, D)
        for j in 1:D
            X[:, j] = chebPts(2π*rand(N), a[j], b[j])
        end

        return X
    end

    return a'  .+ (b - a)' .* rand(N, D)
end

function initializeSol(D::Int, a::Vector{Float64}, b::Vector{Float64})
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

function getBest(fitness::Vector{Float64}, searchType::Symbol = :minimize)
    if searchType == :minimize
        best_X = indmin(fitness) # minimization.
        best = fitness[best_X] 
    else
        best_X = indmax(fitness) # maximization.
        best = fitness[best_X] 
    end

    return best_X, best
end

function Selection(fOld::Individual_, fNew::Individual_, searchType::Symbol)
    gOld = fOld.g .> 0
    hOld = fOld.h .!= 0

    gNew = fNew.g .> 0
    hNew = fNew.h .!= 0

    if sum(gNew) + sum(hNew) < sum(gOld) + sum(hOld)
        return true
    end

    if searchType == :minimize
        return fNew.f < fOld.f
    end
    
    return fNew.f > fOld.f
end