abstract type AbstractAlgorithm end

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


struct Problem <: AbstractProblem
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
