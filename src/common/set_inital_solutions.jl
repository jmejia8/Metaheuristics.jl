include("gen_initial_state.jl")

function set_inital_solutions!(algo::AbstractAlgorithm, solution::AbstractSolution)
    status = algo.status
    if !isnothing(status.best_sol) && !isempty(status.population)
        push!(status.population, solution)
    else
        algo.status = State(solution, [solution])
    end
    algo
end


function set_inital_solutions!(algo::AbstractAlgorithm, x::AbstractVector, fx)
    set_inital_solutions!(algo, create_child(x, fx))
end

function set_inital_solutions!(algo::AbstractAlgorithm, x::AbstractVector, f::Function)
    set_inital_solutions!(algo, x, f(x))
end


function set_inital_solutions!(algo::AbstractAlgorithm, X::AbstractMatrix, fX::AbstractVector)
    n = size(X, 1)
    m = length(fX)

    if n != m
        @warn "$(n) decision vectors provided but $(m) objective values."
        n = min(m, n)
        println("Taking ", n, " as the number of initial solutions.")
    end

    # nothing to do due to it is necessary the objective value
    n == 0 && (return algo)

    # TODO: this part can be parallelized
    for i in 1:n
        set_inital_solutions!(algo, X[i,:], fX[i])
    end

    # TODO check population size provided in algo.parameters.N

    algo
end

function set_inital_solutions!(algo::AbstractAlgorithm, X::AbstractMatrix, f::Function)
    set_inital_solutions!(algo, X, [f(X[i,:]) for i in 1:size(X,1)])
end

