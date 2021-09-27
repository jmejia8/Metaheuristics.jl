mutable struct xf_indiv <: AbstractSolution # Single Objective
    x::Vector{Float64}
    f::Float64
end

mutable struct xfgh_indiv <: AbstractSolution # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
    is_feasible::Bool

end

function xfgh_indiv(
    x::Vector{Float64},
    f::Float64,
    g::Vector{Float64},
    h::Vector{Float64};
    sum_violations = 0.0,
    ε = 0.0
)
    if sum_violations <= 0.0
        sum_violations = violationsSum(g, h; ε=ε)
    end

    xfgh_indiv(x, f, g, h, sum_violations, sum_violations == 0.0)
end

mutable struct xFgh_indiv <: AbstractSolution# Single Objective Constraied
    x::Vector{Float64}
    f::Vector{Float64}
    g::Vector{Float64}
    h::Vector{Float64}
    rank::Int
    crowding::Float64
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
    is_feasible::Bool
end

function xFgh_indiv(
    x::Vector{Float64},
    f::Vector{Float64},
    g::Vector{Float64},
    h::Vector{Float64};
    rank = 0,
    crowding = 0.0,
    sum_violations = 0.0,
    ε = 0.0
)

    if sum_violations <= 0
        sum_violations = violationsSum(g, h;ε=ε)
    end
    xFgh_indiv(x, f, g, h, Int(rank), crowding, sum_violations, sum_violations == 0.0)
end



const Solution = Union{xf_indiv, xfgh_indiv, xFgh_indiv}
const Population = Array{Solution, 1}


############################################################
# Generate solutions depending on the objective function
# output
#############################################################

function generateChild(x::Vector{Float64}, fResult::Float64;ε=0.0)
    return xf_indiv(x, fResult)
end

function generateChild(x::Vector{Float64},
        fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}};
        ε=0.0
    )
    f, g, h = fResult
    return xfgh_indiv(x, f, g, h; ε=ε)
end

function generateChild(x::Vector{Float64},
        fResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}};
        ε = 0.0
    )
    f, g, h = fResult
    return xFgh_indiv(x, f, g, h;ε=ε)
end



function generate_population(N::Int, problem;ε=0.0)
    a = view(problem.bounds, 1, :)'
    b = view(problem.bounds, 2, :)'
    D = length(a)

    X = a .+ (b - a) .* rand(N, D)

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
create_child(x, fx) = generateChild(x, fx)

function create_solution(x, problem::Problem; ε=0.0)
    problem.f_calls += 1
    return generateChild(x, problem.f(x), ε=ε)
end



# getters for the above structures

"""
    get_position(solution)

Get the position vector.
"""
get_position(solution::Solution) = solution.x

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
fval(solution::Solution) = solution.f
sum_violations(solution::xfgh_indiv) = solution.sum_violations
sum_violations(solution::xFgh_indiv) = solution.sum_violations


function fvals(population::AbstractArray{xFgh_indiv})
    if isempty(population)
        return zeros(0, 0)
    end
    
    M = length(fval(population[1]))
    [ fval(sol)[i] for sol in population, i in 1:M]
end

fvals(population::AbstractArray) = fval.(population)



"""
    pareto_front(state::State)
Returns the non-dominated solutions in state.population.
"""
pareto_front(st::State) = pareto_front(st.population)


"""
    pareto_front(population::Array)
Returns non-dominated solutions.
"""
pareto_front(population::Array) = fvals(get_non_dominated_solutions(population))

