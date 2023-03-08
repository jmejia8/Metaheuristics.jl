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
  iteration: 500
    minimum: 8.63808e-08
  minimizer: [0.0002658098418993323, -1.140808975532608e-5, -0.00012488307670533095]
    f calls: 15000
 total time: 0.1556 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], CGSA(N = 80, chaosIndex = 1))
+=========== RESULT ==========+
  iteration: 500
    minimum: 1.0153e-09
  minimizer: [-8.8507563788141e-6, -1.3050111801923072e-5, 2.7688577445980026e-5]
    f calls: 40000
 total time: 1.0323 s
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

    # random initialization for agents.
    status = gen_initial_state(problem,parameters,information,options,status)
    N = parameters.N
    options.f_calls_limit = options.f_calls_limit == 0 ? options.iterations * N : options.f_calls_limit

    # Velocity
    parameters.V = isempty(parameters.V) ? zeros(N,D) : parameters.V
    # Postions
    parameters.X = isempty(parameters.X) ? positions(status) : parameters.X
    # function values
    if isempty(parameters.fitness)
        parameters.fitness = fvals(status.population)
    end
    

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
    end

    status.population = create_solutions(X, problem)

    parameters.X = X
    parameters.fitness = fvals(status.population)
    parameters.V = V

    #Evaluation of agents. 
    currentBest = get_best(status.population)

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

