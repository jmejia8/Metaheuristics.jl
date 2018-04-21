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

function initializePop(func::Function, N::Int, D::Int, a::Vector{Float64}, b::Vector{Float64}, initType::Symbol=:uniform)
    X = initializePop(N, D, a, b, initType)
    
    # infers datatype
    x = X[1,:]
    child = generateChild(x, func(x))
    individual = typeof(child)

    # population array
    population = Array{individual, 1}([])

    # first individual
    push!(population, child)

    for i in 2:N
        x = X[i,:]
        child = generateChild(x, func(x))
        push!(population, child)
    end

    return population
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
    v = fobj(X[1,:])
    soltype = typeof(v)

    f = Array{soltype}([])
    push!(f, v)
    for i = 2:N
        push!(f, fobj(X[i,:]))
    end

    return f
end

function getfValues(f::Vector{Float64})
    return f
end

function getfValues(fv::Array{Tuple{Float64,Array{Float64,1}}})
    f = zeros(length(fv))

    for i = 1:length(fv)
        f[i] = fv[i][1]
    end
    return f
end

function getfValues(fv::Array{ Tuple{Float64,Array{Float64,1},Array{Float64,1}} })
    f = zeros(length(fv))

    for i = 1:length(fv)
        f[i] = fv[i][1]
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

# COP functions
function countViolations(g::Vector, h::Vector)
    sum_g = 0
    sum_h = 0

    for i = 1:length(g)
        if g[i] > 0
            sum_g += g[i]  end
    end

    for i = 1:length(h)
        if h[i] != 0.0
            sum_h += abs(h[i])  end
    end

    return sum_g + sum_h
end
