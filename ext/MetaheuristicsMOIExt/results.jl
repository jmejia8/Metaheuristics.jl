function MOI.get(model::Optimizer, ::MOI.SolveTimeSec)
    isnothing(model.results) && return 0.0
    return model.results.overall_time
end

function converged(res)
    res.termination_status_code isa Metaheuristics.ConvergenceTermination
end

_best_is_feasible(res) = !isempty(res.population) &&
    Metaheuristics.is_feasible(res.best_sol)

function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    isnothing(model.infeasible_message) || return MOI.INFEASIBLE
    isnothing(model.results) && return MOI.OPTIMIZE_NOT_CALLED

    term = model.results.termination_status_code

    if term isa _SETTLED
        # settled on an infeasible point: no feasible point was found
        return _found_feasible(model) ? MOI.LOCALLY_SOLVED : MOI.LOCALLY_INFEASIBLE
    elseif term isa Metaheuristics.IterationLimit
        return MOI.ITERATION_LIMIT
    elseif term isa Metaheuristics.EvaluationsLimit
        return MOI.OTHER_LIMIT
    elseif term isa Metaheuristics.TimeLimit
        return MOI.TIME_LIMIT
    else
        # For other termination reasons (e.g., OtherLimit, UnknownStopReason)
        isempty(model.results.population) && return MOI.OTHER_LIMIT
        return _found_feasible(model) ? MOI.LOCALLY_SOLVED : MOI.LOCALLY_INFEASIBLE
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    isnothing(model.infeasible_message) || return model.infeasible_message
    isnothing(model.results) && return "optimize! not called"
    msg = Metaheuristics.termination_status_message(model.results)
    isnothing(model.pareto) && return msg
    return "$msg Found $(length(model.pareto)) solution(s)."
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    isnothing(model.results) && return 0
    isnothing(model.pareto) || return length(model.pareto)
    1
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end

    if !isnothing(model.pareto)
        point = model.pareto[attr.result_index]
        return point.feasible ? MOI.FEASIBLE_POINT : MOI.INFEASIBLE_POINT
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
    isnothing(model.pareto) || return copy(model.pareto[attr.result_index].y)
    val = minimum(model.results)
    if model.sense == MOI.MAX_SENSE
        val = -val
    end
    val
end

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    isnothing(model.pareto) || return model.pareto[attr.result_index].x[vi.value]
    solution = Metaheuristics.vec_to_dict(minimizer(model.results), model.search_space)

    x = _dict_to_values(solution, model.search_space, model.var_dist)
    x[vi.value]
end

# Multi-objective: the ideal point of the reported solutions (best value of each
# objective). As with MultiObjectiveAlgorithms, it is not a proven bound here.
function MOI.get(model::Optimizer, attr::MOI.ObjectiveBound)
    if isnothing(model.pareto) || isempty(model.pareto)
        throw(MOI.GetAttributeNotAllowed(
            attr,
            "ObjectiveBound is only available after a multi-objective solve.",
        ))
    end
    best = model.sense == MOI.MAX_SENSE ? max : min
    return reduce((a, b) -> best.(a, b), (p.y for p in model.pareto))
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
