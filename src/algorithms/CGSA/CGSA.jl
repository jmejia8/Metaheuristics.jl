# GSA code v0.1.
# Coded by Jesús Mejía. 
# Based on MATLAB code of Esmat Rashedi, 2010. 
# Adopted from 
# https://la.mathworks.com/matlcgsaentral/fileexchange/61116-gsa-+-chaotic-gravitational-constant
# " E. Rashedi, H. Nezamabadi-pour and S. Saryazdi,
# “GSA: A Gravitational Search Algorithm”, Information sciences, vol. 179,
# no. 13, pp. 2232-2248, 2009."
# ------------------------------------------
# Mirjalili, Seyedali, and Amir H. Gandomi. 
# "Chaotic gravitational constants for the gravitational search algorithm." 
# Applied Soft Computing 53 (2017): 407-419.
# ------------------------------------------

include("physics.jl")
include("chaos.jl")

mutable struct CGSA <: AbstractParameters
    N::Int    
    chValueInitial::Real   
    chaosIndex::Real   
    ElitistCheck::Int    
    Rpower::Int    
    Rnorm::Int    
    wMax::Real   
    wMin::Real   
    X::Matrix{Float64}
    V::Matrix{Float64}
    fitness::Vector{Float64}
end

"""

    CGSA(;
        N::Int    = 30,
        chValueInitial::Real   = 20,
        chaosIndex::Real   = 9,
        ElitistCheck::Int    = 1,
        Rpower::Int    = 1,
        Rnorm::Int    = 2,
        wMax::Real   = chValueInitial,
        wMin::Real   = 1e-10,
        information = Information(),
        options = Options()
    )

CGSA is an extension of the GSA algorithm but with Chaotic gravitational constants for
the gravitational search algorithm.

Ref. Chaotic gravitational constants for the gravitational search algorithm.
Applied Soft Computing 53 (2017): 407-419.

Parameters:


- N: Population size
- chValueInitial: Initial value for the chaos value
- chaosIndex: Integer 1 ≤ chaosIndex ≤ 10 is the function that model the chaos
- Rpower: power related to the distance norm(x)^Rpower
- Rnorm: is the value as in norm(x, Rnorm)

# Example


```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], CGSA())
+=========== RESULT ==========+
| Iter.: 499
| f(x) = 0.000235956
| solution.x = [0.0028549782101697785, -0.0031385153631797724, 0.014763299731686608]
| f calls: 15000
| Total time: 0.1003 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], CGSA(N = 80, chaosIndex = 1))
+=========== RESULT ==========+
| Iter.: 499
| f(x) = 0.000102054
| solution.x = [0.00559987302269564, 0.00017535321765604905, 0.008406213942044265]
| f calls: 40000
| Total time: 0.5461 s
+============================+
```


"""
function CGSA(;
        N::Int    = 30,
        chValueInitial::Real   = 20,
        chaosIndex::Real   = 8,
        ElitistCheck::Int    = 1,
        Rpower::Int    = 1,
        Rnorm::Int    = 2,
        wMax::Real   = chValueInitial,
        wMin::Real   = 1e-10,
        X::Matrix{Float64} = zeros(0,0),
        V::Matrix{Float64} = zeros(0,0),
        fitness::Vector{Float64} = zeros(0),
        information = Information(),
        options = Options()
    )

    parameters = CGSA(N, chValueInitial, chaosIndex, ElitistCheck, Rpower, Rnorm, wMax, wMin, X, V, fitness)


    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end

function initialize!(
        status,
        parameters::CGSA,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
       )


    Rnorm = parameters.Rnorm
    N = parameters.N
    D = size(problem.bounds, 2)
    fobj = problem.f


    # bounds vectors
    low, up = problem.bounds[1,:], problem.bounds[2,:]

    max_it = 500
    options.iterations = options.iterations == 0 ? max_it : options.iterations
    options.f_calls_limit = options.f_calls_limit == 0 ? options.iterations * N : options.f_calls_limit

    # random initialization for agents.
    P = generate_population(N, problem)

    # Current best
    theBest = get_best(P)


    status = State(theBest, P)
    status.f_calls = N

    # Velocity
    parameters.V = isempty(parameters.V) ? zeros(N,D) : parameters.V
    # Postions
    parameters.X = isempty(parameters.X) ? positions(status) : parameters.X
    # function values
    parameters.fitness = isempty(parameters.fitness) ? fvals(status) : parameters.fitness

    return status

end

function update_state!(
        status,
        parameters::CGSA,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
        )


    wMax = parameters.wMax
    N = parameters.N
    wMin = parameters.wMin
    max_it = options.iterations
    iteration = status.iteration
    searchType = :minimize
    chaosIndex = parameters.chaosIndex	
    Rnorm = parameters.Rnorm
    Rpower = parameters.Rpower
    ElitistCheck = parameters.ElitistCheck

    X = parameters.X
    V = parameters.V
    fitness = parameters.fitness

    P = status.population
    theBest = status.best_sol

    # iteration
    chValue = wMax-iteration*((wMax-wMin)/max_it)

    #Calculation of M. eq.14-20
    M = massCalculation(fitness,searchType)

    #Calculation of Gravitational constant. eq.13.
    G = Gconstant(iteration, max_it)

    if 1 <= chaosIndex <= 10
        G += chaos(chaosIndex,iteration,max_it,chValue)
    end


    #Calculation of accelaration in gravitational field. eq.7-10,21.
    a = Gfield(M,X,G,Rnorm,Rpower,ElitistCheck,iteration,max_it)

    #Agent movement. eq.11-12
    X, V = move(X,a,V)


    #
    # Checking allowable range. 
    # X = correctPop(X, low, up)
    for i = 1:N
        x = reset_to_violated_bounds!(X[i,:], problem.bounds)
        X[i,:] = x
        P[i] = create_solution(x, problem)
        fitness[i] = P[i].f
    end
    status.f_calls = problem.f_calls

    parameters.X = X
    parameters.fitness = fitness
    parameters.V = V
    status.population = P

    #Evaluation of agents. 
    currentBest = get_best(P)

    # fix this
    if is_better(currentBest, theBest)
        status.best_sol = currentBest
    end

    stop_criteria!(status, parameters, problem, information, options)

end


function final_stage!(
        status,
        parameters::CGSA,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()
end

