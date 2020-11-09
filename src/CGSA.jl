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



mutable struct CGSA
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
        chaosIndex::Real   = 9,
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
        initialize! = initialize_cgsa!,
        update_state! = update_state_cgsa!,
        is_better = is_better_eca,
        stop_criteria = stop_check,
        final_stage! = final_stage_cgsa!,
        information = information,
        options = options,
    )
end

function initialize_cgsa!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
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
    P = initializePop(fobj, N, D, low, up)
    status.population = P
    status.f_calls = N

    # Current best
    theBest = getBest(P, :minimize)
    status.best_sol = theBest

    # Velocity
    parameters.V = isempty(parameters.V) ? zeros(N,D) : parameters.V
    # Postions
    parameters.X = isempty(parameters.X) ? positions(status) : parameters.X
    # function values
    parameters.fitness = isempty(parameters.fitness) ? fvals(status) : parameters.fitness

end

function update_state_cgsa!(
        problem,
        engine,
        parameters,
        status,
        information,
        options,
        iteration,
        )


    wMax = parameters.wMax
    N = parameters.N
    wMin = parameters.wMax
    max_it = options.iterations
    searchType = :minimize
    chaosIndex = parameters.chaosIndex	
    Rnorm = parameters.Rnorm
    Rpower = parameters.Rpower
    ElitistCheck = parameters.ElitistCheck
    low, up = problem.bounds[1,:], problem.bounds[2,:]

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

    # Checking allowable range. 
    X = correctPop(X, low, up)
    for i = 1:N
        x = X[i,:]
        P[i] = generateChild(x, problem.f(x))
        fitness[i] = P[i].f
    end
    status.f_calls += N

    parameters.X = X
    parameters.fitness = fitness
    parameters.V = V
    status.population = P

    #Evaluation of agents. 
    currentBest = getBest(P, :minimize)

    # fix this
    if engine.is_better(currentBest, theBest)
        status.best_sol = currentBest
    end

    status.stop = engine.stop_criteria(status, information, options)

end


function final_stage_cgsa!(status, information, options)
    status.final_time = time()
end





function chaos(index,curr_iter,max_iter,Value)
    x = zeros(max_iter + 1)
    x[1]=0.7

    G = zeros(max_iter)
    if index == 1
        # Chebyshev map
        for i=1:max_iter
            x[i+1]=cos(i*acos(x[i]))
            G[i]=((x[i]+1)*Value)/2
        end
    elseif index == 2
        # Circle map
        a=0.5
        b=0.2
        for i=1:max_iter
            x[i+1]=mod(x[i]+b-(a/(2*pi))*sin(2*pi*x[i]),1)
            G[i]=x[i]*Value
        end
    elseif index == 3
        # Gauss/mouse map
        for i=1:max_iter
            if x[i]==0
                x[i+1]=0
            else
                x[i+1]=mod(1/x[i],1)
            end
            G[i]=x[i]*Value
        end
    elseif index == 4
        # Iterative map
        a=0.7
        for i=1:max_iter
            x[i+1]=sin((a*pi)/x[i])
            G[i]=((x[i]+1)*Value)/2
        end
    elseif index == 5
        # Logistic map
        a=4
        for i=1:max_iter
            x[i+1]=a*x[i]*(1-x[i])
            G[i]=x[i]*Value
        end
    elseif index == 6
        # Piecewise map
        P=0.4
        for i=1:max_iter
            if x[i]>=0 && x[i]<P
                x[i+1]=x[i]/P
            end
            if x[i]>=P && x[i]<0.5
                x[i+1]=(x[i]-P)/(0.5-P)
            end
            if x[i]>=0.5 && x[i]<1-P
                x[i+1]=(1-P-x[i])/(0.5-P)
            end
            if x[i]>=1-P && x[i]<1
                x[i+1]=(1-x[i])/P
            end    
            G[i]=x[i]*Value
        end

    elseif index == 7
        # Sine map
        for i=1:max_iter
            x[i+1] = sin(pi*x[i])
            G[i]=(x[i])*Value
        end
    elseif index == 8
        # Singer map 
        u=1.07
        for i=1:max_iter
            x[i+1] = u*(7.86*x[i]-23.31*(x[i]^2)+28.75*(x[i]^3)-13.302875*(x[i]^4))
            G[i]=(x[i])*Value
        end
    elseif index == 9
        # Sinusoidal map
        for i=1:max_iter
            x[i+1] = 2.3*x[i]^2*sin(pi*x[i])
            G[i]=(x[i])*Value
        end

    elseif index == 10
        # Tent map
        x[1]=0.6
        for i=1:max_iter
            if x[i]<0.7
                x[i+1]=x[i]/0.7
            end
            if x[i]>=0.7
                x[i+1]=(10/3)*(1-x[i])
            end
            G[i]=(x[i])*Value
        end

    end
    return G[curr_iter]

end

