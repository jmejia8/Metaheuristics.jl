abstract type AbstractSolution end
abstract type AbstractAlgorithm end

mutable struct xf_indiv <: AbstractSolution # Single Objective
    x::Vector{Float64}
    f::Float64
end

mutable struct xfg_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
end

mutable struct xfgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|

end

function xfgh_indiv(
    x::Vector{Float64},
    f::Float64,
    g::Vector{Float64},
    h::Vector{Float64};
    sum_violations = 0
)
    if sum_violations < 0
        sum_violations = violationsSum(g, h)
    end

    xfgh_indiv(x, f, g, h, sum_violations)
end

mutable struct xFgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Vector{Float64}
    g::Vector{Float64}
    h::Vector{Float64}
    rank::Int
    crowding::Float64
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
end

function xFgh_indiv(
    x::Vector{Float64},
    f::Vector{Float64},
    g::Vector{Float64},
    h::Vector{Float64};
    rank = 0,
    crowding = 0.0,
    sum_violations = 0.0
)

    if sum_violations <= 0
        sum_violations = violationsSum(g, h)
    end
    xFgh_indiv(x, f, g, h, Int(rank), crowding, sum_violations)
end

mutable struct State
    best_sol::Any
    population::Array

    f_calls::Int
    g_calls::Int
    h_calls::Int

    iteration::Int
    success_rate::Float64
    convergence::Array{State}

    start_time::Float64
    final_time::Float64
    stop::Bool

end

function State(
    best_sol,
    population;

    # upper level parameters
    f_calls = 0,
    g_calls = 0,
    h_calls = 0,
    iteration = 0,
    success_rate = 0,
    convergence = State[],
    start_time = 0.0,
    final_time = 0.0,
    stop = false,
)

    State(#
        best_sol,
        Array(population),

        # upper level parameters
        promote(f_calls, g_calls, h_calls, iteration)...,
        Real(success_rate),
        State[],
        start_time,
        final_time,
        stop,
    )

end

mutable struct Engine
    initialize!::Function
    update_state!::Function
    is_better::Function
    stop_criteria::Function
    final_stage!::Function
end

function Engine(;
    initialize!::Function = _1(kwargs...) = nothing,
    update_state!::Function = _2(kwargs...) = nothing,
    is_better::Function = _4(kwargs...) = false,
    stop_criteria::Function = _5(kwargs...) = nothing,
    final_stage!::Function = _6(kwargs...) = nothing,
)

    Engine(initialize!, update_stat, is_better, stop_criteria, final_stage!)
end




#####################################################
#
#         STRUCTURES FOR THE OPTIONS
#
#####################################################

mutable struct Options
    x_tol::Float64
    f_tol::Float64
    g_tol::Float64
    h_tol::Float64
    f_calls_limit::Float64
    g_calls_limit::Float64
    h_calls_limit::Float64

    iterations::Int
    store_convergence::Bool
    show_results::Bool
    debug::Bool
    search_type::Symbol
end

function Options(;
    x_tol::Real = 0.0,
    f_tol::Real = 0.0,
    g_tol::Real = 0.0,
    h_tol::Real = 0.0,
    f_calls_limit::Real = 0,
    g_calls_limit::Real = 0,
    h_calls_limit::Real = 0,
    iterations::Int = 0,
    store_convergence::Bool = false,
    show_results::Bool = true,
    debug::Bool = false,
    search_type::Symbol = :minimize,
)


    Options(
        promote(x_tol, f_tol, g_tol, h_tol)...,
        promote(f_calls_limit, g_calls_limit, h_calls_limit)...,
        promote(iterations)...,

        # Results options
        promote(store_convergence, show_results, debug)...,
        Symbol(search_type),
    )

end



struct Information
    f_optimum::Float64
    x_optimum::Array{Float64}
end

function Information(;#
    f_optimum = NaN,
    x_optimum::Array{Float64} = Float64[],
)

    Information(Float64(f_optimum), x_optimum)

end


mutable struct Algorithm <: AbstractAlgorithm
    parameters::Any
    status::State
    information::Information
    options::Options
    engine::Engine
end

function Algorithm(
    parameters;
    initial_state::State = State(nothing, []),
    initialize!::Function = _1(kwargs...) = nothing,
    update_state!::Function = _2(kwargs...) = nothing,
    # is_better(a, b)  = true if x is better that y
    is_better::Function = _5(kwargs...) = false,
    stop_criteria::Function = stop_check,
    final_stage!::Function = _4(kwargs...) = nothing,
    information::Information = Information(),
    options::Options = Options(),
)


    engine = Engine(
        initialize!,
        update_state!,
        is_better,
        stop_criteria,
        final_stage!,
    )

    Algorithm(parameters, initial_state, information, options, engine)

end


struct Problem
    f::Function
    bounds::Array{Float64,2}
    g::Array{Function}
    h::Array{Function}
    type::Symbol
end

function Problem(f::Function, bounds::Array, g = Function[], h = Function[])

    type::Symbol = :constrained

    if length(g) == 0 && length(h) == 0
        type = :unconstrained
    else
        type = :constrained
    end

    Problem(f, bounds, g, h, type)
end
