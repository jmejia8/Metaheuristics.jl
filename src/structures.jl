abstract type AbstractSolution end

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
end

mutable struct xFgh_indiv # Single Objective Constraied
	x::Vector{Float64}
	f::Vector{Float64}
	g::Vector{Float64}
	h::Vector{Float64}
end

mutable struct State{T<:AbstractSolution}
    best_sol::T
    population::Array{T,1}

    f_calls::Int
    g_calls::Int
    h_calls::Int

    iteration::Int
    success_rate::Float64
    convergence::Array{State}

end

function State(
        best_sol,
        population;

        # upper level parameters
        f_calls = 0,
        g_calls = 0,
        h_calls = 0,
        
        iteration= 0,

        success_rate= 0,
        convergence = State[],
    )

    State{xf_indiv}(#
        best_sol,
        Array(population),
        
        # upper level parameters
        promote(
            f_calls,
            g_calls,
            h_calls,
            iteration)...,
            Real(success_rate),
            State[])
    
end

mutable struct Engine
    initialize!::Function
    update_state!::Function
    is_better::Function
    stop_criteria::Function
    final_stage!::Function
end

function Engine(;initialize!::Function = _1(kwargs...) = nothing,
                   update_state!::Function = _2(kwargs...) = nothing,
                       is_better::Function = _4(kwargs...) = false, 
                   stop_criteria::Function = _5(kwargs...) = nothing,
                    final_stage!::Function = _6(kwargs...) = nothing)
    
    Engine(initialize!,update_stat,is_better,stop_criteria,final_stage!)
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
    
    iterations::Int = 1000,
    store_convergence::Bool = false,
    show_results::Bool = true,
    debug::Bool = false,
    search_type::Symbol=:minimize)

    
    Options(
        promote(x_tol, f_tol, g_tol, h_tol)...,
        promote(f_calls_limit, g_calls_limit, h_calls_limit)...,
        
        promote(iterations)...,

        # Results options
        promote(store_convergence,show_results, debug)...,
        Symbol(search_type)
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


mutable struct Algorithm
    parameters
    status::State
    information::Information
    options::Options
    engine::Engine
end

function Algorithm(   parameters;
                   initial_state::State    = State(nothing, []),
                      initialize!::Function = _1(kwargs...) = nothing,
                   update_state!::Function = _2(kwargs...) = nothing,
                       # is_better(a, b)  = true if x is better that y 
                       is_better::Function = _5(kwargs...) = false,
                   stop_criteria::Function = stop_check,
                    final_stage!::Function = _4(kwargs...) = nothing,
                     information::Information = Information(),
                         options::Options  = Options())
    

    engine = Engine(initialize!,
                update_state!,
                is_better,
                stop_criteria,
                final_stage!)

    Algorithm(  parameters,
                initial_state,
                information,
                options,
                engine)

end


struct Problem
    f::Function
    bounds::Matrix{Float64}
    g::Function
    h::Function
    type::Symbol
end

function Problem(F::Function,
                f::Function,
                bounds::Array,
                g::Function = Function[],
                h::Function = Function[])

    type::Symbol = :constrained

    if length(g) == 0 && length(h) == 0
        type = :unconstrained
    else
        type = :constrained
    end
    
    Problem(f, bounds, g, h, type)
end


