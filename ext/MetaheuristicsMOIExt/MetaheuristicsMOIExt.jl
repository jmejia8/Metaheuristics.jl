module MetaheuristicsMOIExt

using Metaheuristics
using SearchSpaces
import MathOptInterface as MOI

function __init__()
    @static if isdefined(Base, :get_extension)
        setglobal!(Metaheuristics, :Optimizer, Optimizer)
    end
    println("Extension was loaded!")
end

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Problem data.
    variables::MOI.Utilities.VariablesContainer{T}
    starting_values::Vector{Union{Nothing, T}}
    nlp_model::Union{MOI.Nonlinear.Model, Nothing}
    # legacy nonlinear interface (`@NLobjective`, `@NLconstraint`)
    nlp_block::Union{MOI.NLPBlockData, Nothing}
    sense::MOI.OptimizationSense

    # Parameters.
    method::Union{Metaheuristics.AbstractAlgorithm, Nothing}
    silent::Bool
    options::Dict{Symbol, Any}

    # Solution attributes.
    results::Union{Nothing, Metaheuristics.State}
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
        MOI.FEASIBILITY_SENSE,
        nothing,
        false,
        Dict{Symbol, Any}(),
        nothing,
        nothing,
        Dict(),
        nothing,
        nothing,
        x -> x,
    )
end

Optimizer() = Optimizer{Float64}()

# Equality constraints h(x) = 0 are handled as |h(x)| - tol <= 0, where tol is
# the `h_tol` of the method options, or this value when `h_tol` is zero.
const DEFAULT_H_TOL = 1.0e-6

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

function MOI.set(model::Optimizer, ::MOI.NLPBlock, data::Union{Nothing, MOI.NLPBlockData})
    model.nlp_block = data
    return
end

MOI.get(model::Optimizer, ::MOI.NLPBlock) = model.nlp_block

function MOI.supports(::Optimizer, ::Union{MOI.ObjectiveSense, MOI.ObjectiveFunction})
    return true
end
MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.supports(::Optimizer, p::MOI.RawOptimizerAttribute)
    return p.name == "method" || hasfield(Metaheuristics.Options, Symbol(p.name))
end

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex})
    # TODO: starting solution
    return false
end

# TODO: add ZeroOne constraints, and integer one
const BOUNDS{T} = Union{
    MOI.LessThan{T},
    MOI.GreaterThan{T},
    MOI.EqualTo{T},
    MOI.Interval{T},
}
const _SETS{T} = Union{MOI.GreaterThan{T}, MOI.LessThan{T}, MOI.EqualTo{T}}

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
    # TODO add compat with MOI.VectorNonlinearFunction
    nl = convert(MOI.ScalarNonlinearFunction, func)

    if isnothing(model.nlp_model)
        model.nlp_model = MOI.Nonlinear.Model()
    end
    MOI.Nonlinear.set_objective(model.nlp_model, nl)
    nothing
end

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    model.silent = value
    nothing
end

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

const TIME_LIMIT = "time_limit"
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
    MOI.supports(model, p) || throw(MOI.UnsupportedAttribute(p))
    key = Symbol(p.name)
    # convert now, so a wrong value fails here instead of inside `optimize!`
    model.options[key] = convert(fieldtype(Metaheuristics.Options, key), value)
    return
end

# Copy the user options (`RawOptimizerAttribute`, `TimeLimitSec`, `Silent`)
# into the options of the method used for solving.
function _apply_options!(opts::Metaheuristics.Options, model::Optimizer)
    for (key, value) in model.options
        setfield!(opts, key, value)
    end
    # a new seed needs a new generator, unless the user also gave one
    if haskey(model.options, :seed) && !haskey(model.options, :rng)
        opts.rng = Metaheuristics.default_rng_mh(opts.seed)
    end
    # Silent = false keeps the verbosity chosen by the user
    model.silent && (opts.verbose = false)
    return opts
end

function MOI.get(model::Optimizer, p::MOI.RawOptimizerAttribute)
    if p.name == "method"
        return model.method
    end
    key = Symbol(p.name)
    if haskey(model.options, key)
        return model.options[key]
    end
    error("RawOptimizerAttribute with name $(p.name) is not set.")
end

function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    isnothing(model.results) && return 0.0
    return model.results.overall_time
end

function MOI.empty!(model::Optimizer)
    MOI.empty!(model.variables)
    empty!(model.starting_values)
    model.nlp_model = nothing
    model.nlp_block = nothing
    model.sense = MOI.FEASIBILITY_SENSE
    model.results = nothing
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

function starting_value(optimizer::Optimizer{T}, i) where {T}
    if optimizer.starting_values[i] !== nothing
        return optimizer.starting_values[i]
    else
        v = optimizer.variables
        return min(max(zero(T), v.lower[i]), v.upper[i])
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real, Nothing},
)
    MOI.throw_if_not_valid(model, vi)
    model.starting_values[vi.value] = value
    return
end

function _dict_to_search_spaces(dict)

    s = Pair[]
    var_dist = Dict()

    if dict[:binary][:count] > 0
        push!(s, :binary => BitArraySpace(dict[:binary][:count]))
        var_dist[:binary] = dict[:binary][:vars]
    end

    if !isempty(dict[:int][:lb])
        push!(s, :int => BoxConstrainedSpace(dict[:int][:lb], dict[:int][:ub]))
        var_dist[:int] = dict[:int][:vars]
    end

    if !isempty(dict[:bounds][:lb])
        push!(s, :bounds => BoxConstrainedSpace(dict[:bounds][:lb], dict[:bounds][:ub]))
        var_dist[:bounds] = dict[:bounds][:vars]
    end

    if !isempty(dict[:non_rigid][:int][:lb])
        lb = dict[:non_rigid][:int][:lb]
        ub = dict[:non_rigid][:int][:ub]
        push!(s, :int_non_rigid => BoxConstrainedSpace(lb, ub, rigid = false))
        var_dist[:int_non_rigid] = dict[:non_rigid][:int][:vars]
    end

    if !isempty(dict[:non_rigid][:bounds][:lb])
        lb = dict[:non_rigid][:bounds][:lb]
        ub = dict[:non_rigid][:bounds][:ub]
        push!(s, :bounds_non_rigid => BoxConstrainedSpace(lb, ub, rigid = false))
        var_dist[:bounds_non_rigid] = dict[:non_rigid][:bounds][:vars]
    end

    MixedSpace(s...), var_dist
end

# Integer bounds at or beyond this value are treated as infinite in the search
# space (an `Int` box cannot hold them, and floats lose integer precision).
const INT_BOUND_LIMIT = 2.0^52

# Return a message if the variable bounds admit no point, otherwise `nothing`.
function _check_bounds(model)
    lb = model.variables.lower
    ub = model.variables.upper
    for vi in MOI.get(model.variables, MOI.ListOfVariableIndices())
        a, b = lb[vi.value], ub[vi.value]
        k = get(model.variable_type, vi, :bounds)
        if k === :binary
            a, b = max(a, 0), min(b, 1)
        end
        if k !== :bounds
            a, b = ceil(a), floor(b)
        end
        if a > b
            kind = k === :bounds ? "" : k === :int ? "integer " : "binary "
            return "Infeasible bounds: $(kind)variable $(vi.value) has no value " *
                "in [$(lb[vi.value]), $(ub[vi.value])]."
        end
    end
    return nothing
end

function _get_search_spaces(model)

    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    lb = model.variables.lower |> copy
    ub = model.variables.upper |> copy

    spaces = Dict(
        :binary => Dict(:count => 0, :vars => []),
        :non_rigid => Dict(
            :int => Dict(:lb => Int[], :ub => Int[], :vars => []),
            :bounds => Dict(:lb => Float64[], :ub => Float64[], :vars => []),
        ),
        :int => Dict(:lb => Int[], :ub => Int[], :vars => []),
        :bounds => Dict(:lb => Float64[], :ub => Float64[], :vars => []),
    )

    for (a, b, i) in zip(lb, ub, vars)
        # Get variable type, default to :bounds if not set
        k = get(model.variable_type, i, :bounds)

        if k === :binary
            spaces[k][:count] += 1
            push!(spaces[k][:vars], i)
            continue
        end

        if k === :int
            # Too large bounds cannot be stored as Int: handle them as infinite
            # here (they are still enforced as penalty constraints).
            abs(a) < INT_BOUND_LIMIT || (a = -Inf)
            abs(b) < INT_BOUND_LIMIT || (b = Inf)
        end

        if isfinite(a) && isfinite(b)
            if k === :int
                # fractional bounds: tightest integer box
                a, b = ceil(Int, a), floor(Int, b)
            end
            push!(spaces[k][:lb], a)
            push!(spaces[k][:ub], b)
            push!(spaces[k][:vars], i)
            continue
        end

        if isinf(a) && isinf(b)
            # TODO:  study the performance on different initial bound values on
            # unbounded search spaces
            a, b = -10, 10
        elseif isinf(a)
            a = b - 2max(abs(b), 10)
        else
            b = a + 2max(abs(a), 10)
        end
        if k === :int
            a, b = ceil(Int, a), floor(Int, b)
        end
        push!(spaces[:non_rigid][k][:lb], a)
        push!(spaces[:non_rigid][k][:ub], b)
        push!(spaces[:non_rigid][k][:vars], i)
    end
    spaces
end

function _dict_to_values(solution, search_space, var_dist)
    # TODO improve this part
    x = zeros(Real, Metaheuristics.getdim(search_space))
    for (k, vs) in var_dist
        for (i, v) in enumerate(vs)
            x[v.value] = solution[k][i]
        end
    end
    x
end

# Nonlinear functions come from two places: the new interface
# (`ScalarNonlinearFunction`, stored in `model.nlp_model`) and the legacy
# `NLPBlock` (`@NLobjective`, `@NLconstraint`). Both evaluators take `x` in MOI
# variable order. Only zero-order evaluations are needed here.
function _nonlinear_sources(model::Optimizer, vars)
    sources = Tuple{Any, Vector{MOI.NLPBoundsPair}}[]
    objective = nothing

    if !isnothing(model.nlp_model)
        backend = MOI.Nonlinear.SparseReverseMode()
        evaluator = MOI.Nonlinear.Evaluator(model.nlp_model, backend, vars)
        MOI.initialize(evaluator, Symbol[])
        bounds = MOI.NLPBlockData(evaluator).constraint_bounds
        isempty(bounds) || push!(sources, (evaluator, bounds))
        isnothing(model.nlp_model.objective) || (objective = evaluator)
    end

    if !isnothing(model.nlp_block)
        block = model.nlp_block
        MOI.initialize(block.evaluator, Symbol[])
        isempty(block.constraint_bounds) ||
            push!(sources, (block.evaluator, block.constraint_bounds))
        # as in other MOI solvers, an NLPBlock objective takes precedence
        block.has_objective && (objective = block.evaluator)
    end

    # FEASIBILITY_SENSE ignores the objective; no objective means f(x) = 0
    model.sense == MOI.FEASIBILITY_SENSE && (objective = nothing)
    return objective, sources
end

function MOI.optimize!(model::Optimizer{T}) where {T}
    model.results = nothing
    model.infeasible_message = _check_bounds(model)
    if !isnothing(model.infeasible_message)
        !model.silent && @warn model.infeasible_message
        return
    end

    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    objective_scale = model.sense == MOI.MAX_SENSE ? -one(T) : one(T)
    obj_evaluator, nl_sources = _nonlinear_sources(model, vars)

    # constraint bounds in the order of the concatenated constraint values
    nl_bounds = reduce(vcat, (b for (_, b) in nl_sources); init = MOI.NLPBoundsPair[])
    nl_constrained = !isempty(nl_bounds)
    lc = [b.lower for b in nl_bounds]
    uc = [b.upper for b in nl_bounds]
    eq_idx = findall(lc .== uc)
    up_idx = findall(i -> lc[i] != uc[i] && isfinite(uc[i]), eachindex(uc))
    lo_idx = findall(i -> lc[i] != uc[i] && isfinite(lc[i]), eachindex(lc))

    !model.silent && @info "Load search space"
    search_space_dict = _get_search_spaces(model)
    # get search space and the variable mapping JuMP <--> Metaheuristics
    search_space, var_dist = _dict_to_search_spaces(search_space_dict)
    model.search_space = search_space
    model.var_dist = var_dist

    # Variable bounds that the search space does not enforce by itself, so they
    # must be handled as constraints g(x) <= 0: the finite bound of semi-bounded
    # variables (non-rigid spaces) and bounds on binary variables (BitArraySpace
    # ignores them). Bounds of variables in rigid boxes are always satisfied.
    lb = model.variables.lower
    ub = model.variables.upper
    bound_constraints = Tuple{Int, Symbol, T}[] # (variable index, :lower/:upper, bound)
    for key in (:int_non_rigid, :bounds_non_rigid, :binary)
        for vi in get(var_dist, key, ())
            a, b = lb[vi.value], ub[vi.value]
            is_bin = key === :binary
            # [0, 1] bounds are implied for binary variables
            isfinite(a) &&
                (!is_bin || a > 0) &&
                push!(bound_constraints, (vi.value, :lower, a))
            isfinite(b) &&
                (!is_bin || b < 1) &&
                push!(bound_constraints, (vi.value, :upper, b))
        end
    end
    has_constraints = nl_constrained || !isempty(bound_constraints)

    !model.silent && @info "Load method"

    # Solve with a fresh copy: `Algorithm` is mutable and `optimize` continues
    # from `method.status` (population and fitness from a previous solve, maybe
    # of another objective or dimension), and it also mutates parameters and
    # termination criteria. `model.method` is left exactly as the user set it.
    method = if isnothing(model.method)
        # TODO add compat with multiobjective optimization
        options = Options(verbose = !model.silent, iterations = 10_000, f_calls_limit = Inf)
        MixedInteger(Metaheuristics.ECA(adaptive = true); options)
    else
        deepcopy(model.method)
    end
    method.status = Metaheuristics.State(nothing, [])
    _apply_options!(method.options, model)

    # Equalities are handled as |h(x)| - h_tol <= 0 inside g, since not every
    # algorithm uses `options.h_tol` when it decides feasibility (ECA does not).
    h_tol = method.options.h_tol > 0 ? method.options.h_tol : DEFAULT_H_TOL

    C = zeros(length(nl_bounds))
    function constraints(x)
        GV = T[
            side === :lower ? bound - x[i] : x[i] - bound
            for (i, side, bound) in bound_constraints
        ]
        nl_constrained || return GV, zeros(1)

        offset = 0
        for (evaluator, bounds) in nl_sources
            c = zeros(length(bounds))
            MOI.eval_constraint(evaluator, c, x)
            C[(offset + 1):(offset + length(c))] = c
            offset += length(c)
        end

        G = [
            C[up_idx] - uc[up_idx]                # g(x) <= u
            lc[lo_idx] - C[lo_idx]                # g(x) >= l
            GV                                    # variable bounds
            abs.(C[eq_idx] - uc[eq_idx]) .- h_tol # h(x) = b
        ]
        return G, zeros(1)
    end

    # TODO: multi-objective
    function f(solution)
        x = _dict_to_values(solution, search_space, var_dist)
        fx = if isnothing(obj_evaluator)
            zero(T)
        else
            objective_scale * MOI.eval_objective(obj_evaluator, x)
        end

        has_constraints || return fx
        gx, hx = constraints(x)
        fx, gx, hx
    end

    model.f = f
    !model.silent && @info "Optimize"

    # TODO initial solutions
    #initial_x = starting_value.(model, eachindex(model.starting_values))
    model.results = optimize(f, search_space, method)

    return
end

function converged(res)
    res.termination_status_code isa Metaheuristics.ConvergenceTermination
end

_best_is_feasible(res) = !isempty(res.population) &&
    Metaheuristics.is_feasible(res.best_sol)

# Stops meaning "the search settled" (not "the budget ran out").
const _SETTLED = Union{
    Metaheuristics.ConvergenceTermination,
    Metaheuristics.ObjectiveVarianceLimit,
    Metaheuristics.ObjectiveDifferenceLimit,
    Metaheuristics.AccuracyLimit,
}

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    isnothing(model.infeasible_message) || return MOI.INFEASIBLE
    isnothing(model.results) && return MOI.OPTIMIZE_NOT_CALLED

    term = model.results.termination_status_code

    if term isa _SETTLED
        # settled on an infeasible point: no feasible point was found
        return _best_is_feasible(model.results) ?
            MOI.LOCALLY_SOLVED :
            MOI.LOCALLY_INFEASIBLE
    elseif term isa Metaheuristics.IterationLimit
        return MOI.ITERATION_LIMIT
    elseif term isa Metaheuristics.EvaluationsLimit
        return MOI.OTHER_LIMIT
    elseif term isa Metaheuristics.TimeLimit
        return MOI.TIME_LIMIT
    else
        # For other termination reasons (e.g., OtherLimit, UnknownStopReason)
        isempty(model.results.population) && return MOI.OTHER_LIMIT
        return _best_is_feasible(model.results) ?
            MOI.LOCALLY_SOLVED :
            MOI.LOCALLY_INFEASIBLE
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    isnothing(model.infeasible_message) || return model.infeasible_message
    isnothing(model.results) && return "optimize! not called"
    return Metaheuristics.termination_status_message(model.results)
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    isnothing(model.results) && return 0
    # TODO: check for multi-objective case
    # return length(model.results.population)
    1
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    # Check if the best solution is feasible
    if !isempty(model.results.population)
        best = model.results.best_sol
        if Metaheuristics.is_feasible(best)
            return MOI.FEASIBLE_POINT
        else
            return MOI.INFEASIBLE_POINT
        end
    end

    return MOI.UNKNOWN_RESULT_STATUS
end
MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    val = minimum(model.results)
    if model.sense == MOI.MAX_SENSE
        val = -val
    end
    val
end

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    solution = Metaheuristics.vec_to_dict(minimizer(model.results), model.search_space)

    x = _dict_to_values(solution, model.search_space, model.var_dist)
    x[vi.value]
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{
        MOI.VariableIndex,
        <:Union{BOUNDS{T}, MOI.ZeroOne, MOI.Integer},
    },
) where {T}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    # The minimizer is stored in the MixedSpace order; map back to MOI order.
    vi = MOI.VariableIndex(ci.value)
    return MOI.get(model, MOI.VariablePrimal(attr.result_index), vi)
end
end # module
