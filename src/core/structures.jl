abstract type AbstractAlgorithm end


mutable struct Algorithm{T} <: AbstractAlgorithm
    parameters::T
    status::State
    information::Information
    options::Options
end

function Algorithm(
    parameters;
    initial_state::State = State(nothing, []),
    information::Information = Information(),
    options::Options = Options(),
)

    Algorithm(parameters, initial_state, information, options)

end

function Base.show(io::IO, alg::Algorithm)
    Base.show(io, alg.parameters)
end

termination_status_message(alg::Algorithm) = termination_status_message(alg.status)
should_stop(algorithm::AbstractAlgorithm) = algorithm.status.stop
function get_result(algorithm::AbstractAlgorithm)
    if isnothing(algorithm.status.best_sol)
        error("First optimize a function: `optimize!(f, bounds, method)`")
    end
    
    algorithm.status
end

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
    Problem(f, BitArrays(size(bounds,2)); kargs...)
end

function evaluate(x, problem::Problem)
    problem.f_calls += 1
    return problem.f(x)
end

