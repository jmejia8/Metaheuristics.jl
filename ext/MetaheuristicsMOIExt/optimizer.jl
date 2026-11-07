# One reported solution of a multi-objective solve: `x` in MOI variable order and
# `y` the objective vector with its true sign (as in MultiObjectiveAlgorithms).
struct SolutionPoint{T}
    x::Vector{T}
    y::Vector{T}
    feasible::Bool
end

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Problem data.
    variables::MOI.Utilities.VariablesContainer{T}
    starting_values::Vector{Union{Nothing, T}}
    nlp_model::Union{MOI.Nonlinear.Model, Nothing}
    # legacy nonlinear interface (`@NLobjective`, `@NLconstraint`)
    nlp_block::Union{MOI.NLPBlockData, Nothing}
    # vector objective (multi-objective): one constraint row per objective
    obj_model::Union{MOI.Nonlinear.Model, Nothing}
    n_objectives::Int
    sense::MOI.OptimizationSense

    # Parameters.
    method::Union{Metaheuristics.AbstractAlgorithm, Nothing}
    silent::Bool
    options::Dict{Symbol, Any}
    # number of initial solutions built from the start values (`nothing`: default)
    start_seeds::Union{Nothing, Int}

    # Solution attributes.
    results::Union{Nothing, Metaheuristics.State}
    # non-dominated solutions of a multi-objective solve (`nothing` otherwise)
    pareto::Union{Nothing, Vector{SolutionPoint{T}}}
    # set when infeasibility is detected before solving (e.g. empty bounds)
    infeasible_message::Union{Nothing, String}
    variable_type::Dict
    search_space
    var_dist
    f
end

function Optimizer{T}() where {T}
    Optimizer{T}(
        MOI.Utilities.VariablesContainer{T}(),
        Union{Nothing, T}[],
        nothing,
        nothing,
        nothing,
        0,
        MOI.FEASIBILITY_SENSE,
        nothing,
        false,
        Dict{Symbol, Any}(),
        nothing,
        nothing,
        nothing,
        nothing,
        Dict(),
        nothing,
        nothing,
        x -> x,
    )
end

Optimizer() = Optimizer{Float64}()

function MOI.supports(::Optimizer, ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction})
    return true
end
MOI.supports(::Optimizer, ::MOI.Silent) = true

function MOI.supports(::Optimizer, p::MOI.RawOptimizerAttribute)
    return p.name == "method" ||
        p.name == START_SEEDS ||
        hasfield(Metaheuristics.Options, Symbol(p.name))
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.VariableIndex},
    ::Type{<:BOUNDS{T}},
) where {T}
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{MOI.ZeroOne},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{MOI.Integer},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:_SETS{T}},
) where {T}
    return true
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Metaheuristics"

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.sense = sense
    return
end

function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F}
    nl = convert(MOI.ScalarNonlinearFunction, func)

    if isnothing(model.nlp_model)
        model.nlp_model = MOI.Nonlinear.Model()
    end
    MOI.Nonlinear.set_objective(model.nlp_model, nl)
    # a scalar objective replaces a previous vector one
    model.obj_model = nothing
    model.n_objectives = 0
    nothing
end

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    model.silent = value
    nothing
end

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true

function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawOptimizerAttribute(TIME_LIMIT), Float64(value))
end

function MOI.set(model::Optimizer, attr::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options, Symbol(TIME_LIMIT))
end

function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, Symbol(TIME_LIMIT), nothing)
end

MOI.Utilities.map_indices(::Function, opt::Metaheuristics.AbstractAlgorithm) = opt

function MOI.set(model::Optimizer, p::MOI.RawOptimizerAttribute, value)
    if p.name == "method"
        model.method = value
        return
    end
    if p.name == START_SEEDS
        if !(isnothing(value) || (value isa Integer && value >= 1))
            throw(ArgumentError("\"$START_SEEDS\" must be a positive integer or nothing."))
        end
        model.start_seeds = value
        return
    end
    MOI.supports(model, p) || throw(MOI.UnsupportedAttribute(p))
    key = Symbol(p.name)
    # convert now, so a wrong value fails here instead of inside `optimize!`
    model.options[key] = convert(fieldtype(Metaheuristics.Options, key), value)
    return
end

function MOI.get(model::Optimizer, p::MOI.RawOptimizerAttribute)
    if p.name == "method"
        return model.method
    end
    if p.name == START_SEEDS && !isnothing(model.start_seeds)
        return model.start_seeds
    end
    key = Symbol(p.name)
    if haskey(model.options, key)
        return model.options[key]
    end
    error("RawOptimizerAttribute with name $(p.name) is not set.")
end

function MOI.empty!(model::Optimizer)
    MOI.empty!(model.variables)
    empty!(model.starting_values)
    model.nlp_model = nothing
    model.nlp_block = nothing
    model.obj_model = nothing
    model.n_objectives = 0
    model.sense = MOI.FEASIBILITY_SENSE
    model.results = nothing
    model.pareto = nothing
    model.infeasible_message = nothing
    empty!(model.variable_type)
    model.search_space = nothing
    model.var_dist = nothing
    return
end

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.variables) &&
        isempty(model.starting_values) &&
        isnothing(model.nlp_model) &&
        isnothing(model.nlp_block) &&
        isnothing(model.obj_model) &&
        model.sense == MOI.FEASIBILITY_SENSE
end

function MOI.add_variable(model::Optimizer{T}) where {T}
    push!(model.starting_values, nothing)
    vi = MOI.add_variable(model.variables)
    # Set default type for variables without explicit type constraints
    # This will be overridden if a ZeroOne, Integer, or other constraint is added
    model.variable_type[vi] = :bounds
    return vi
end

function MOI.is_valid(
    model::Optimizer,
    index::Union{MOI.VariableIndex, MOI.ConstraintIndex},
)
    return MOI.is_valid(model.variables, index)
end

function MOI.add_constraint(model::Optimizer, vi::MOI.VariableIndex, set::MOI.ZeroOne)
    model.variable_type[vi] = :binary
    # TODO: revise this part to get g(x) <= 0
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    vi::MOI.VariableIndex,
    set::BOUNDS{T},
) where {T}
    # bounds do not change the variable type (it may already be :int or :binary)
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.add_constraint(model::Optimizer, vi::MOI.VariableIndex, set::MOI.Integer)
    # Mark variable as integer (binary is stronger, keep it whatever the order)
    if get(model.variable_type, vi, :bounds) !== :binary
        model.variable_type[vi] = :int
    end
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    f::MOI.ScalarNonlinearFunction,
    s::_SETS{T},
) where {T}
    if isnothing(model.nlp_model)
        model.nlp_model = MOI.Nonlinear.Model()
    end
    index = MOI.Nonlinear.add_constraint(model.nlp_model, f, s)
    MOI.ConstraintIndex{typeof(f), typeof(s)}(index.value)
end
