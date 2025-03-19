mutable struct CSO <: AbstractParameters
    N::Int
    φ::Float64
    V_loser::Matrix{Float64}
end


"""
    CSO(;
        N = 400,
        φ = 0.1,
        information = Information(),
        options = Options()
    )

A competitive swarm optimizer for large scale optimization.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], CSO())
Optimization Result
===================
  Iteration:       149
  Minimum:         3.85394e-18
  Minimizer:       [6.08039e-10, 1.85731e-9, -1.86081e-10]
  Function calls:  30000
  Total time:      0.0382 s
  Stop reason:     Maximum objective function calls exceeded.


julia> optimize(f, [-1 -1 -1; 1 1 1.0], CSO(N = 100, φ = 0.4))
Optimization Result
===================
  Iteration:       476
  Minimum:         3.76679e-51
  Minimizer:       [-4.68574e-26, 3.31086e-26, 2.17945e-26]
  Function calls:  23850
  Total time:      0.0065 s
  Stop reason:     Due to Convergence Termination criterion.
```

"""
function CSO(; N::Int = 400, φ = 0.1, kargs...)
    parameters = CSO(N, φ, zeros(0,0))
    Algorithm( parameters; kargs...)
end

function update_state!(
        status,
        parameters::CSO,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    N = parameters.N
    D = getdim(problem)
    φ = parameters.φ
    V_loser  = parameters.V_loser
    population = status.population

    # Determine the losers and winners
    rank    = randperm(options.rng, N)
    loser   = rank[1:N ÷ 2]
    winner  = rank[(end-N÷2)+1:end]
    mask = [is_better(l, w) for (l, w) in zip(population[loser], population[winner])]
    loser[mask], winner[mask] = winner[mask], loser[mask]
    # Update the losers by learning from the winners
    X_loser  = population[loser] |> positions
    X_winner = population[winner]|> positions
    R1 = rand(options.rng, N ÷ 2, D)
    R2 = rand(options.rng, N ÷ 2, D)
    R3 = rand(options.rng, N ÷ 2, D)
    # velocity of loser particles
    x_mean = mean(positions(population), dims=1)
    parameters.V_loser = R1.*V_loser + R2.*(X_winner-X_loser) + φ *R3.*(x_mean .- X_loser);
    X_loser = X_loser + parameters.V_loser;

    # repair particle positions
    for x in eachrow(X_loser)
        reset_to_violated_bounds!(x, problem.search_space)
    end
    

    population[loser] = create_solutions(X_loser, problem,ε=options.h_tol)
    status.best_sol = get_best(population)
end

function initialize!(
    status,
    parameters::CSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
       )
    D = getdim(problem)


    if parameters.N == 0
        parameters.N = clamp(10 * D, 10, 500)
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
        @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N÷2) + 1
    end

    parameters.V_loser = zeros(parameters.N÷2, getdim(problem))
    gen_initial_state(problem,parameters,information,options,status)
end

function final_stage!(
    status,
    parameters::CSO,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )

end
