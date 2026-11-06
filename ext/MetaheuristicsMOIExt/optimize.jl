"""
    _NonlinearConstraints

Nonlinear constraints `l <= c(x) <= u` of every nonlinear source, laid out in the
order of their concatenated values.

# Fields
- `sources`: `(evaluator, bounds)` pairs, as returned by `_nonlinear_sources`.
- `lower`, `upper`: bounds `l` and `u` of the concatenated constraint values.
- `eq_idx`: indices of equalities (`l == u`).
- `up_idx`: indices of finite upper bounds that are not equalities.
- `lo_idx`: indices of finite lower bounds that are not equalities.
"""
struct _NonlinearConstraints{S}
    sources::S
    lower::Vector{Float64}
    upper::Vector{Float64}
    eq_idx::Vector{Int}
    up_idx::Vector{Int}
    lo_idx::Vector{Int}
end

"""
    _NonlinearConstraints(sources) -> _NonlinearConstraints

Split the bounds of the nonlinear `sources` into equality, upper and lower parts.
"""
function _NonlinearConstraints(sources)
    # constraint bounds in the order of the concatenated constraint values
    bounds = reduce(vcat, (b for (_, b) in sources); init = MOI.NLPBoundsPair[])
    lc = [b.lower for b in bounds]
    uc = [b.upper for b in bounds]
    eq_idx = findall(lc .== uc)
    up_idx = findall(i -> lc[i] != uc[i] && isfinite(uc[i]), eachindex(uc))
    lo_idx = findall(i -> lc[i] != uc[i] && isfinite(lc[i]), eachindex(lc))
    return _NonlinearConstraints(sources, lc, uc, eq_idx, up_idx, lo_idx)
end

Base.length(nl::_NonlinearConstraints) = length(nl.lower)
Base.isempty(nl::_NonlinearConstraints) = isempty(nl.lower)

"""
    _reset_results!(model::Optimizer) -> Optimizer

Discard the results of a previous solve.
"""
function _reset_results!(model::Optimizer)
    model.results = nothing
    model.pareto = nothing
    return model
end

"""
    _has_feasible_bounds!(model::Optimizer) -> Bool

Store in `model.infeasible_message` why the variable bounds are infeasible (or
`nothing`), warn about it unless silent, and return whether they are feasible.
"""
function _has_feasible_bounds!(model::Optimizer)
    model.infeasible_message = _check_bounds(model)
    isnothing(model.infeasible_message) && return true
    !model.silent && @warn model.infeasible_message
    return false
end

"""
    _objective_scale(model::Optimizer{T}) -> T

`-1` for maximization (Metaheuristics minimizes), `1` otherwise.
"""
_objective_scale(model::Optimizer{T}) where {T} = model.sense == MOI.MAX_SENSE ?
    -one(T) :
    one(T)

"""
    _load_search_space!(model::Optimizer) -> (search_space, var_dist)

Build the search space and the variable mapping JuMP <--> Metaheuristics, and
store both in `model`.
"""
function _load_search_space!(model::Optimizer)
    search_space_dict = _get_search_spaces(model)
    search_space, var_dist = _dict_to_search_spaces(search_space_dict)
    model.search_space = search_space
    model.var_dist = var_dist
    return search_space, var_dist
end

"""
    _bound_constraints(model::Optimizer{T}, var_dist) -> Vector{Tuple{Int, Symbol, T}}

Variable bounds that the search space does not enforce by itself, so they must
be handled as constraints g(x) <= 0, as `(variable index, :lower/:upper, bound)`.
These are the finite bound of semi-bounded variables (non-rigid spaces) and
bounds on binary variables (`BitArraySpace` ignores them). Bounds of variables
in rigid boxes are always satisfied.
"""
function _bound_constraints(model::Optimizer{T}, var_dist) where {T}
    lb = model.variables.lower
    ub = model.variables.upper
    bound_constraints = Tuple{Int, Symbol, T}[]
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
    return bound_constraints
end

"""
    _default_method(multiobjective::Bool, silent::Bool) -> MixedInteger

Method used when the user sets none: NSGA2 for vector objectives, adaptive ECA
otherwise.
"""
function _default_method(multiobjective::Bool, silent::Bool)
    if multiobjective
        # NSGA2 picks its own budget (500 iterations) from these options
        options = Options(verbose = !silent, iterations = 5_000, f_calls_limit = Inf)
        return MixedInteger(Metaheuristics.NSGA2(); options)
    end
    options = Options(verbose = !silent, iterations = 10_000, f_calls_limit = Inf)
    return MixedInteger(Metaheuristics.ECA(adaptive = true); options)
end

"""
    _solver_method(model::Optimizer, multiobjective::Bool) -> AbstractAlgorithm

Fresh method to solve `model` with: a copy of `model.method` (or the default
method), with an empty state and the user options applied.

`Algorithm` is mutable and `optimize` continues from `method.status`
(population and fitness from a previous solve, maybe of another objective or
dimension), and it also mutates parameters and termination criteria, so
`model.method` is left exactly as the user set it.
"""
function _solver_method(model::Optimizer, multiobjective::Bool)
    method = if isnothing(model.method)
        _default_method(multiobjective, model.silent)
    else
        _check_method(model.method, multiobjective)
        deepcopy(model.method)
    end
    method.status = Metaheuristics.State(nothing, [])
    _apply_options!(method.options, model)
    return method
end

"""
    _equality_tolerance(options::Options) -> Real

Tolerance of the equality constraints: `options.h_tol`, or `DEFAULT_H_TOL` when
it is not positive.
"""
_equality_tolerance(options) = options.h_tol > 0 ? options.h_tol : DEFAULT_H_TOL

"""
    _bound_violations(::Type{T}, bound_constraints, x) -> Vector{T}

Values g(x) <= 0 of the variable `bound_constraints` at `x`.
"""
_bound_violations(::Type{T}, bound_constraints, x) where {T} = T[
    side === :lower ? bound - x[i] : x[i] - bound
    for (i, side, bound) in bound_constraints
]

"""
    _eval_nonlinear!(C, nl::_NonlinearConstraints, x) -> C

Write the concatenated values of the nonlinear constraints at `x` into `C`.
"""
function _eval_nonlinear!(C, nl::_NonlinearConstraints, x)
    offset = 0
    for (evaluator, bounds) in nl.sources
        c = zeros(length(bounds))
        MOI.eval_constraint(evaluator, c, x)
        C[(offset + 1):(offset + length(c))] = c
        offset += length(c)
    end
    return C
end

"""
    _constraint_function(::Type{T}, nl, bound_constraints, h_tol) -> Union{Function, Nothing}

Function `x -> (g, h)` with every constraint as g(x) <= 0 (and a dummy `h`), or
`nothing` when the problem is unconstrained.

Equalities are handled as |h(x)| - h_tol <= 0 inside g, since not every
algorithm uses `options.h_tol` when it decides feasibility (ECA does not).
"""
function _constraint_function(
    ::Type{T},
    nl::_NonlinearConstraints,
    bound_constraints,
    h_tol,
) where {T}
    nl_constrained = !isempty(nl)
    nl_constrained || !isempty(bound_constraints) || return nothing

    lc, uc = nl.lower, nl.upper
    eq_idx, up_idx, lo_idx = nl.eq_idx, nl.up_idx, nl.lo_idx
    C = zeros(length(nl)) # reused between evaluations
    return function constraints(x)
        GV = _bound_violations(T, bound_constraints, x)
        nl_constrained || return GV, zeros(1)

        _eval_nonlinear!(C, nl, x)
        G = [
            C[up_idx] - uc[up_idx]                # g(x) <= u
            lc[lo_idx] - C[lo_idx]                # g(x) >= l
            GV                                    # variable bounds
            abs.(C[eq_idx] - uc[eq_idx]) .- h_tol # h(x) = b
        ]
        return G, zeros(1)
    end
end

"""
    _single_objective_function(obj_evaluator, objective_scale::T, constraints, search_space, var_dist) -> Function

Objective `solution -> f` (or `(f, g, h)` when `constraints` is not `nothing`)
for single-objective algorithms. No objective evaluator means f(x) = 0.
"""
function _single_objective_function(
    obj_evaluator,
    objective_scale::T,
    constraints,
    search_space,
    var_dist,
) where {T}
    return function f_single(solution)
        x = _dict_to_values(solution, search_space, var_dist)
        fx = if isnothing(obj_evaluator)
            zero(T)
        else
            objective_scale * MOI.eval_objective(obj_evaluator, x)
        end

        isnothing(constraints) && return fx
        gx, hx = constraints(x)
        fx, gx, hx
    end
end

"""
    _multi_objective_function(model::Optimizer{T}, vars, objective_scale::T, constraints, search_space, var_dist) -> Function

Objective `solution -> (f, g, h)` for multi-objective algorithms, evaluating the
vector objective stored in `model.obj_model`.
"""
function _multi_objective_function(
    model::Optimizer{T},
    vars,
    objective_scale::T,
    constraints,
    search_space,
    var_dist,
) where {T}
    mo_evaluator = MOI.Nonlinear.Evaluator(
        model.obj_model,
        MOI.Nonlinear.SparseReverseMode(),
        vars,
    )
    MOI.initialize(mo_evaluator, Symbol[])
    F = zeros(T, model.n_objectives)
    return function f_multi(solution)
        x = _dict_to_values(solution, search_space, var_dist)
        MOI.eval_constraint(mo_evaluator, F, x)
        fx = objective_scale .* F # new vector: F is reused
        # multi-objective algorithms expect (f, g, h)
        isnothing(constraints) && return fx, zeros(1), zeros(1)
        gx, hx = constraints(x)
        fx, gx, hx
    end
end

function MOI.optimize!(model::Optimizer{T}) where {T}
    _reset_results!(model)
    _has_feasible_bounds!(model) || return

    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    objective_scale = _objective_scale(model)
    obj_evaluator, nl_sources = _nonlinear_sources(model, vars)
    nl = _NonlinearConstraints(nl_sources)

    !model.silent && @info "Load search space"
    search_space, var_dist = _load_search_space!(model)
    bound_constraints = _bound_constraints(model, var_dist)

    !model.silent && @info "Load method"
    multiobjective = _is_multiobjective(model)
    method = _solver_method(model, multiobjective)

    h_tol = _equality_tolerance(method.options)
    constraints = _constraint_function(T, nl, bound_constraints, h_tol)
    f = if multiobjective
        _multi_objective_function(
            model,
            vars,
            objective_scale,
            constraints,
            search_space,
            var_dist,
        )
    else
        _single_objective_function(
            obj_evaluator,
            objective_scale,
            constraints,
            search_space,
            var_dist,
        )
    end

    model.f = f
    # after resetting `method.status`, which would discard them
    _set_start_solutions!(method, model, f, search_space, var_dist)

    !model.silent && @info "Optimize"
    model.results = optimize(f, search_space, method)
    if multiobjective
        model.pareto = _pareto_points(model, model.results, objective_scale)
    end

    return
end
