
# Base.@kwdef struct Termination
#     criteria_all::Vector = Any[]
#     criteria_any::Vector = Any[]
# end


mutable struct Problem{S} <: AbstractProblem
    f::Function
    search_space::S
    f_calls::Int
    parallel_evaluation::Bool
end

function Base.getproperty(obj::Problem, sym::Symbol)
    if sym === :bounds
        @warn "`bounds` property is deprecated. Use `search_space` instead"
        se = obj.search_space
        return Array([se.lb se.ub]')
    else # fallback to getfield
        return getfield(obj, sym)
    end
end

function getdim(problem::Problem)
    return SearchSpaces.getdim(problem.search_space)
end


function Problem(f::Function, search_space::AbstractSearchSpace; parallel_evaluation=false)
    Problem(f, search_space, 0, parallel_evaluation)
end

function Problem(f::Function, bounds::AbstractMatrix{T}; kargs...) where T <: Number
    # old problem definition

    if size(bounds,1) > 2 && size(bounds,2) == 2
        bounds = Array(bounds')
    end

    Problem(f, BoxConstrainedSpace(lb = bounds[1,:], ub = bounds[2,:]); kargs...)
end


function Problem(f::Function, bounds::Array{Bool,2}; kargs...)
    Problem(f, BitArraySpace(size(bounds,2)); kargs...)
end

function evaluate(x, problem::Problem)
    problem.f_calls += 1
    return problem.f(x)
end

