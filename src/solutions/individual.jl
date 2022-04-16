abstract type AbstractUnconstrainedSolution <: AbstractSolution end
abstract type AbstractConstrainedSolution   <: AbstractSolution end
# multi-objective solution are constrained by default
abstract type AbstractMultiObjectiveSolution <: AbstractConstrainedSolution end

mutable struct xf_solution{X} <: AbstractUnconstrainedSolution # Single Objective
    x::X
    f::Float64
end

const xf_indiv = xf_solution{Vector{Float64}}

# Single Objective Constrained
mutable struct xfgh_solution{X} <: AbstractConstrainedSolution 
    x::X
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
    is_feasible::Bool

end

const xfgh_indiv = xfgh_solution{Vector{Float64}}

function xfgh_indiv(
    x,
    f::Real,
    g::Vector,
    h::Vector;
    sum_violations = 0.0,
    ε = 0.0
)
    if sum_violations <= 0.0
        sum_violations = violationsSum(g, h; ε=ε)
    end

    xfgh_solution(x, f, g, h, sum_violations, sum_violations == 0.0)
end

mutable struct xFgh_solution{X} <: AbstractMultiObjectiveSolution
    x::X
    f::Vector{Float64}
    g::Vector{Float64}
    h::Vector{Float64}
    rank::Int
    crowding::Float64
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
    is_feasible::Bool
end

const xFgh_indiv = xFgh_solution{Vector{Float64}}

function xFgh_indiv(
    x,
    f::Vector,
    g::Vector,
    h::Vector;
    rank = 0,
    crowding = 0.0,
    sum_violations = 0.0,
    ε = 0.0
)

    if sum_violations <= 0
        sum_violations = violationsSum(g, h;ε=ε)
    end
    xFgh_solution(x, f, g, h, Int(rank), crowding, sum_violations, sum_violations == 0.0)
end



const Solution = Union{xf_indiv, xfgh_indiv, xFgh_indiv}
const Population = Array{Solution, 1}


############################################################
# Generate solutions depending on the objective function
# output
#############################################################

# single objective
"""
    Metaheuristics.create_child(x, fx)

Constructor for a solution depending on the result of `fx`.

## Example

```
julia> import Metaheuristics

julia> Metaheuristics.create_child(rand(3), 1.0)
| f(x) = 1
| solution.x = [0.2700437125780806, 0.5233263210622989, 0.12871108215859772]

julia> Metaheuristics.create_child(rand(3), (1.0, [2.0, 0.2], [3.0, 0.3]))
| f(x) = 1
| g(x) = [2.0, 0.2]
| h(x) = [3.0, 0.3]
| x = [0.9881102595664819, 0.4816273348099591, 0.7742585077942159]

julia> Metaheuristics.create_child(rand(3), ([-1, -2.0], [2.0, 0.2], [3.0, 0.3]))
| f(x) = [-1.0, -2.0]
| g(x) = [2.0, 0.2]
| h(x) = [3.0, 0.3]
| x = [0.23983577719146854, 0.3611544510766811, 0.7998754930109109]

julia> population = [ Metaheuristics.create_child(rand(2), (randn(2),  randn(2), rand(2))) for i = 1:100  ]
                           F space
          ┌────────────────────────────────────────┐ 
        2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠂⠀⠀⠀⠀⠀⡇⠈⡀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠘⠀⡇⠀⠀⠘⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠂⠀⠀⢀⠠⠐⠀⡇⠄⠁⠀⠀⠀⡀⠀⢁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠂⢈⠀⠈⡇⠀⡐⠃⠀⠄⠄⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠄⢐⠠⠀⡄⠀⠀⡇⠀⠂⠈⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠉⠉⠉⠉⠋⠉⠉⠉⠉⠉⠉⠙⢉⠉⠙⠉⠉⡏⠉⠉⠩⠋⠉⠩⠉⠉⠉⡉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉│ 
   f_2    │⠀⠀⠀⠀⠀⡀⠀⠀⠀⠄⠀⠀⡀⠀⠀⠂⠀⡇⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠐⡇⠠⠀⠀⠀⠈⢀⠄⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠄⠀⡀⠀⠂⡇⠐⠘⠈⠂⠀⠈⡀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠂⠀⠂⠀⠀⡇⠀⠈⢀⠐⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⢁⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       -3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          └────────────────────────────────────────┘ 
          -3                                       4
                             f_1
```

"""
function create_child(x, fResult::Float64;ε=0.0)
    return xf_solution(x, fResult)
end

# constrained single objective
function create_child(x,
        fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}};
        ε=0.0
    )
    f, g, h = fResult
    isempty(g) && isempty(h) && (return xf_solution(x, f))

    return xfgh_indiv(x, f, g, h; ε=ε)
end

# constrained multi-objective
function create_child(x,
        fResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}};
        ε = 0.0
    )
    f, g, h = fResult
    return xFgh_indiv(x, f, g, h;ε=ε)
end

##########################################################3
# parallel evaluation
##########################################################3
# single objective
function create_child(X::AbstractMatrix, fResult::AbstractVector;ε=0.0)
    size(X,1) != length(fResult) && error("Error in parallel evaluation: size(X,2) != length(f(X))")
    return [xf_indiv(X[i,:], fResult[i]) for i in 1:size(X,1)]
end


# unconstrained multi objective
function create_child(X::AbstractMatrix, fResult::AbstractMatrix;ε=0.0)
    size(fResult,2) > 1 && (return create_child(X, (fResult, zeros(0,0), zeros(0,0))))

    size(X,1) != size(fResult,1) && error("Error in parallel evaluation: size(X,2) != length(f(X))") 

    return [xf_indiv(X[i,:], fResult[i]) for i in 1:size(X,1)]
end

# constrained single objective
function create_child(X::AbstractMatrix,
        fResult::Tuple{AbstractVector, T, T};
        ε=0.0
    ) where T <: AbstractMatrix{Float64}

    F, G, H = fResult

    size_ok = size(X,1) == length(F) #&& size(X,1) == size(G,1) && size(X,1) == size(H,1)
    !size_ok && error("Error in parallel evaluation: Size of X, F(X),G(X) or H(X) differs.")

    @assert isempty(G) || !isempty(G) && size(X,1) == size(G,1)
    @assert isempty(H) || !isempty(H) && size(X,1) == size(H,1)

    population = xfgh_indiv[]

    for i in 1:size(X,1)
        g = isempty(G) ? zeros(0) : G[i,:]
        h = isempty(H) ? zeros(0) : H[i,:]
        push!(population, xfgh_indiv(X[i,:], F[i], g, h;ε=ε))
    end
    population
end

# constrained multi-objective
function create_child(X::AbstractMatrix,
        fResult::Tuple{T, T, T};
        ε=0.0
    ) where T <: AbstractMatrix{Float64}

    F, G, H = fResult

    # single-objective
    size(F,2) == 1 && (return create_child(X, (F[:,1], G, H)))

    @assert size(X,1) == size(F,1)
    @assert isempty(G) || !isempty(G) && size(X,1) == size(G,1)
    @assert isempty(H) || !isempty(H) && size(X,1) == size(H,1)
    # !size_ok && error("Error in parallel evaluation: size(X,1) != size(F(X),1).")
    #&& error("Error in parallel evaluation: size(X,1) != size(G(X),1).")
    #&& error("Error in parallel evaluation: size(X,1) != size(H(X),1).")

    population = xFgh_indiv[]

    for i in 1:size(X,1)
        g = isempty(G) ? zeros(0) : G[i,:]
        h = isempty(H) ? zeros(0) : H[i,:]
        push!(population, xFgh_indiv(X[i,:], F[i,:], g, h;ε=ε))
    end
    population
end

##########################################################3

function _gen_X(N, bounds::AbstractMatrix{T}) where T <: AbstractFloat 
    a = view(bounds, 1, :)'
    b = view(bounds, 2, :)'
    D = length(a)
    a .+ (b - a) .* rand(eltype(bounds), N, D)
end

function _gen_X(N, bounds::AbstractMatrix{T}) where T <: Integer
    a = view(bounds, 1, :)'
    b = view(bounds, 2, :)'
    D = length(a)

    X = zeros(eltype(bounds), N, D)
    for i in 1:D
        X[:,i] = rand(a[i]:b[i], N)
    end
    return X
end

function _gen_X(N, bounds::AbstractMatrix{T}) where T <: Bool
    rand(eltype(bounds), N, size(bounds, 2))
end

function generate_population(N::Int, problem;ε=0.0, parallel_evaluation = false)
    X = _gen_X(N, problem.bounds)

    if problem.parallel_evaluation
        return create_solutions(X, problem; ε=ε)
    end
    
    population = [ create_solution(X[i,:], problem; ε=ε) for i in 1:N]

    return population
end


"""
    Metaheuristics.create_child(x, fx)

Constructor for a solution depending on the result of `fx`.

## Example

```
julia> import Metaheuristics

julia> Metaheuristics.create_child(rand(3), 1.0)
| f(x) = 1
| solution.x = [0.2700437125780806, 0.5233263210622989, 0.12871108215859772]

julia> Metaheuristics.create_child(rand(3), (1.0, [2.0, 0.2], [3.0, 0.3]))
| f(x) = 1
| g(x) = [2.0, 0.2]
| h(x) = [3.0, 0.3]
| x = [0.9881102595664819, 0.4816273348099591, 0.7742585077942159]

julia> Metaheuristics.create_child(rand(3), ([-1, -2.0], [2.0, 0.2], [3.0, 0.3]))
| f(x) = [-1.0, -2.0]
| g(x) = [2.0, 0.2]
| h(x) = [3.0, 0.3]
| x = [0.23983577719146854, 0.3611544510766811, 0.7998754930109109]

julia> population = [ Metaheuristics.create_child(rand(2), (randn(2),  randn(2), rand(2))) for i = 1:100  ]
                           F space
          ┌────────────────────────────────────────┐ 
        2 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠂⠀⠀⠀⠀⠀⡇⠈⡀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢀⠘⠀⡇⠀⠀⠘⠀⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠂⠀⠀⢀⠠⠐⠀⡇⠄⠁⠀⠀⠀⡀⠀⢁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠂⢈⠀⠈⡇⠀⡐⠃⠀⠄⠄⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠄⢐⠠⠀⡄⠀⠀⡇⠀⠂⠈⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠉⠉⠉⠉⠋⠉⠉⠉⠉⠉⠉⠙⢉⠉⠙⠉⠉⡏⠉⠉⠩⠋⠉⠩⠉⠉⠉⡉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉⠉│ 
   f_2    │⠀⠀⠀⠀⠀⡀⠀⠀⠀⠄⠀⠀⡀⠀⠀⠂⠀⡇⠀⠀⠀⠐⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠄⠀⠀⠐⡇⠠⠀⠀⠀⠈⢀⠄⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠄⠀⡀⠀⠂⡇⠐⠘⠈⠂⠀⠈⡀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠄⠀⠀⠀⠀⠀⠂⠀⠂⠀⠀⡇⠀⠈⢀⠐⠀⠈⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠁⠀⠀⠀⠀⠀⢁⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠠⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
       -3 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          └────────────────────────────────────────┘ 
          -3                                       4
                             f_1
```

"""
generateChild(x, fx) = create_child(x, fx)

function create_solution(x::AbstractVector, problem::Problem; ε=0.0)
    problem.f_calls += 1
    return create_child(x, problem.f(x), ε=ε)
end


function create_solutions(X::AbstractMatrix, problem::Problem; ε=0.0)
    problem.f_calls += size(X,1)
    if problem.parallel_evaluation
        return create_child(X, problem.f(X), ε=ε)
    end
    [create_child(X[i,:], problem.f(X[i,:]), ε=ε) for i in 1:size(X,1)]
end

# getters for the above structures

"""
    get_position(solution)

Get the position vector.
"""
get_position(solution::AbstractSolution) = solution.x

function positions(population::AbstractArray)
    if isempty(population)
        zeros(0,0)
    end

    D = length(get_position(population[1]))
    [ get_position(sol)[i] for sol in population, i in 1:D]
end
"""
    fval(solution)

Get the objective function value (fitness) of a solution.
"""
fval(solution::AbstractSolution) = solution.f
sum_violations(solution::T) where T <: AbstractConstrainedSolution = solution.sum_violations

"""
    gval(solution)

Get the inequality constraints of a solution.
"""
gval(solution::T) where T <: AbstractConstrainedSolution = solution.g


"""
    hval(solution)

Get the equality constraints of a solution.
"""
hval(solution::T) where T <: AbstractConstrainedSolution = solution.h


function fvals(population::AbstractArray{T}) where T <: AbstractMultiObjectiveSolution
    isempty(population) && (return zeros(0, 0))
    
    M = length(fval(population[1]))
    [ fval(sol)[i] for sol in population, i in 1:M]
end

fvals(population::AbstractArray{T}) where T <: AbstractSolution = fval.(population)


function gvals(population::AbstractArray{T}) where T <: AbstractConstrainedSolution
    if isempty(population)
        return zeros(0, 0)
    end
    
    M = length(gval(population[1]))
    [ gval(sol)[i] for sol in population, i in 1:M]
end


function hvals(population::AbstractArray{T}) where T <: AbstractConstrainedSolution
    if isempty(population)
        return zeros(0, 0)
    end
    
    M = length(hval(population[1]))
    [ hval(sol)[i] for sol in population, i in 1:M]
end


"""
    is_feasible(solution)
Returns `true` if solution is feasible, otherwise returns `false`.
"""
is_feasible(solution::T) where T <: AbstractConstrainedSolution = solution.is_feasible
is_feasible(solution::T) where T <: AbstractSolution = true

"""
    pareto_front(state::State)
Returns the non-dominated solutions in state.population.
"""
pareto_front(st::State) = pareto_front(st.population)


"""
    pareto_front(population::Array)
Returns non-dominated solutions.
"""
function pareto_front(population::AbstractArray{T}) where T <: AbstractMultiObjectiveSolution
    return fvals(get_non_dominated_solutions(population))
end

