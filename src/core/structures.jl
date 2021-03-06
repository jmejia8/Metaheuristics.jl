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

mutable struct Problem <: AbstractProblem
    f::Function
    bounds::Array{Float64,2}
    f_calls::Int
end

function Problem(f::Function, bounds::Array)

    if size(bounds,1) > 2 && size(bounds,2) == 2
        bounds = Array(bounds')
    end

    Problem(f, bounds, 0)
end

function evaluate(x, problem::Problem)
    problem.f_calls += 1
    return problem.f(x)
end


