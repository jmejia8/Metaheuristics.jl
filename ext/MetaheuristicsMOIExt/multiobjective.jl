# Vector objective (`@objective(model, Min, [f1, f2])`): multi-objective problem.
# As in MultiObjectiveAlgorithms, `ObjectiveSense` applies to every component.
# Each component is stored as a free constraint row, so one `eval_constraint`
# call gives the whole objective vector.
function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F <:
    MOI.AbstractVectorFunction}
    obj_model = MOI.Nonlinear.Model()
    rows = MOI.Utilities.scalarize(func)
    for fi in rows
        nl = convert(MOI.ScalarNonlinearFunction, fi)
        MOI.Nonlinear.add_constraint(obj_model, nl, MOI.Interval(-Inf, Inf))
    end
    # a vector objective replaces a previous scalar one
    isnothing(model.nlp_model) || MOI.Nonlinear.set_objective(model.nlp_model, nothing)
    model.obj_model = obj_model
    model.n_objectives = length(rows)
    nothing
end

_base_parameters(p) = p isa MixedInteger ? _base_parameters(p.base) : p
_is_mo_method(method) = _base_parameters(method.parameters) isa _MO_PARAMETERS

# A vector objective is solved as a multi-objective problem, unless the sense is
# FEASIBILITY_SENSE or a legacy `@NLobjective` takes precedence.
function _is_multiobjective(model::Optimizer)
    isnothing(model.obj_model) && return false
    model.sense == MOI.FEASIBILITY_SENSE && return false
    nlp_obj = !isnothing(model.nlp_block) && model.nlp_block.has_objective
    return !nlp_obj
end

function _check_method(method, multiobjective::Bool)
    name = nameof(typeof(_base_parameters(method.parameters)))
    if multiobjective && !_is_mo_method(method)
        throw(ArgumentError(
            "The objective is a vector, but $name is a single-objective " *
                "algorithm. Use a multi-objective one, e.g. MixedInteger(NSGA2()).",
        ))
    elseif !multiobjective && _is_mo_method(method)
        throw(ArgumentError(
            "$name is a multi-objective algorithm, but the objective is scalar. " *
                "Use a vector objective, e.g. @objective(model, Min, [f1, f2]).",
        ))
    end
    return
end

# Non-dominated solutions of the final population, MultiObjectiveAlgorithms
# style: feasible points only (when there is any), duplicated objective vectors
# removed, and sorted lexicographically by objective vector (reversed for Max).
function _pareto_points(model::Optimizer{T}, res, objective_scale) where {T}
    front = Metaheuristics.get_non_dominated_solutions(res.population)
    feasible = filter(Metaheuristics.is_feasible, front)
    isempty(feasible) || (front = feasible)

    points = SolutionPoint{T}[]
    for sol in front
        y = objective_scale .* Metaheuristics.fval(sol)
        any(p -> isapprox(p.y, y; atol = 1.0e-6), points) && continue
        d = Metaheuristics.vec_to_dict(Metaheuristics.get_position(sol), model.search_space)
        x = _dict_to_values(d, model.search_space, model.var_dist)
        push!(points, SolutionPoint{T}(x, y, Metaheuristics.is_feasible(sol)))
    end
    return sort!(points; by = p -> p.y, rev = model.sense == MOI.MAX_SENSE)
end

# Multi-objective: a feasible point was found when any reported point is feasible.
function _found_feasible(model::Optimizer)
    isnothing(model.pareto) || return any(p -> p.feasible, model.pareto)
    return _best_is_feasible(model.results)
end
