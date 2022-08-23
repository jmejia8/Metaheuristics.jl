##################################################################
# Source code taken from: https://github.com/jbytecode/MCCGA
# Credits to corresponding author.
##################################################################
include("utils.jl")

mutable struct MCCGA <: AbstractParameters
    N::Int # population size
    maxsamples::Int
    mutation::Float64
    probvector::Vector{Float64}
    use_local_search::Bool
end

"""
    MCCGA(;N, maxsamples)

### Parameters:

- `N` population size. Default is 100.
- `maxsamples` maximum number of samples. Default is 10000.

### Description 

MCCGA method implements the Machine-coded genetic algorithms for real valued optimization problems. 
The algorithm is based on the concept of a compact genetic algorithm but with a machine-coded representation 
using IEEE-754 floating point encoding standard. In the first stage of the algorithm, `maxsamples` number 
of samples are generated within the range of function domain. This process is required to obtain a vector
of probabilities for each single bit of the IEEE-754 representation. In classical CGAs, the initial vector 
of probabilities is generated using the constant probability of 0.5 whereas in MCCGA, the probability of ith bit having
a value of 1 depends on the function domain. The second step performs a CGA search but with IEEE-754 bits again. Since 
CGA does not use a population of solutions but a single vector of probabilities, the parameter `N` does not really mean 
number of solutions. Instead, it means the amount of mutation in each iteration, e.g. 1/N. 
In the last stage, a local search is performed for fine-tuning. In this implementation, `Optim` package is required.


### References 

- Harik, G. R., Lobo, F. G., & Goldberg, D. E. (1999). The compact genetic algorithm. IEEE transactions on evolutionary computation, 3(4), 287-297.
- Satman, M. H. & Akadal, E. (2020). Machine Coded Compact Genetic Algorithms for Real Parameter Optimization Problems . Alphanumeric Journal , 8 (1) , 43-58 . DOI: 10.17093/alphanumeric.576919
- Mehmet Hakan Satman, Emre Akadal, Makine Kodlu Hibrit Kompakt Genetik Algoritmalar Optimizasyon Yöntemi, Patent, TR, 2022/01, 2018-GE-510239

### Example

```jldoctest

julia> import Optim 

julia> f, bounds, solutions = Metaheuristics.TestProblems.rastrigin();

julia> result = optimize(f, bounds, MCCGA())

+=========== RESULT ==========+
  iteration: 42833
    minimum: 0
  minimizer: [-5.677669786379456e-17, 2.7942451898022582e-39, -3.60925916059986e-33, -6.609510861086017e-34, 2.998586759675011e-32, -1.8825832500775007e-38, -3.0729484147585938e-31, 1.7675578057632446e-38, 5.127944874215823e-16, 1.9623770480857484e-19]
    f calls: 86404
 total time: 3.2839 s
+============================+
```

### Explicit Example 

```jldoctest

julia> using Metaheuristics                                          
julia> using Optim                                                   
                                                              
julia> function f(x::Vector{Float64})::Float64 = (x[1]-pi)^2 + (x[2]-exp(1))^2                                                                                
                                                              
julia> bounds = [ -500.0  -500.0;                                    
                   500.0  500.0]                                      

julia> # Both Metaheuristics and Optim includes optimize() function
julia> result = Metaheuristics.optimize(f, bounds, MCCGA())

+=========== RESULT ==========+
  iteration: 2974
    minimum: 1.80997e-09
  minimizer: [3.1416284249228976, 2.7183048585824263]
    f calls: 6012
 total time: 1.5233 s
stop reason: Other stopping criteria.
+============================+
```
"""
function MCCGA(;
        N = 100,
        maxsamples =10_000,
        mutation = 1 / N,
        use_local_search = true,
        information = Information(),
        options = Options()
    )

    parameters = MCCGA(N, maxsamples, mutation, [], use_local_search)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end

function initialize!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    if options.iterations == 0
        options.iterations = 100_000_000
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end
    
    lower = problem.bounds[1,:]
    upper = problem.bounds[2,:]

    parameters.probvector = initialprobs(lower, upper, maxsamples = parameters.maxsamples)

    # sample a vector to create an initial State
    x = sample(parameters.probvector) |> floats
    initial_sol = create_solution(x, problem)
    return State(initial_sol, [initial_sol for i in 1:parameters.N])

end

function update_state!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    probvector = parameters.probvector
    chsize = length(probvector)

    # parents
    ch1 = sample(probvector)
    ch2 = sample(probvector)

    # evaluate cost function
    sol1 = create_solution(floats(ch1), problem)
    sol2 = create_solution(floats(ch2), problem)

    sol_winner = sol1
    winner = ch1
    loser  = ch2

    # check if ch2 (sol2) is better that ch1 (sol1)
    if is_better(sol2, sol1)
        sol_winner = sol2
        winner = ch2
        loser = ch1
    end

    # save in population the winner (informative only)
    # this is not used in the algorithm
    status.population[1 + status.iteration % parameters.N] = sol_winner

    # save best solution found so far
    is_better(sol_winner, status.best_sol) &&  (status.best_sol = sol_winner)

    for i = 1:chsize
        if winner[i] != loser[i]
            if winner[i] == 1
                probvector[i] += parameters.mutation
            else
                probvector[i] -= parameters.mutation
            end
            if probvector[i] > 1
                probvector[i] = 1
            elseif probvector[i] < 0
                probvector[i] = 0
            end
        end
    end
end

function final_stage!(
        status,
        parameters::MCCGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    status.final_time = time()

    parameters.use_local_search && @warn "The package `Optim` must be installed to perform the local search."
    
    if status.stop
        return
    end
    
end

function stop_criteria!(status, parameters::MCCGA, problem, information, options)
    # check budget limitation
    if status.stop
        return
    end

    mutation = parameters.mutation
    status.stop = all(x -> (x <= mutation) || (x >= 1.0 - mutation), parameters.probvector)
    
    status.stop && (status.termination_status_code = OTHER_LIMIT)
end

