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