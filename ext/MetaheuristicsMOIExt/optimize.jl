function MOI.optimize!(model::Optimizer{T}) where {T}
    model.results = nothing
    model.pareto = nothing
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
    multiobjective = _is_multiobjective(model)
    method = if isnothing(model.method) && multiobjective
        # NSGA2 picks its own budget (500 iterations) from these options
        options = Options(verbose = !model.silent, iterations = 5_000, f_calls_limit = Inf)
        MixedInteger(Metaheuristics.NSGA2(); options)
    elseif isnothing(model.method)
        options = Options(verbose = !model.silent, iterations = 10_000, f_calls_limit = Inf)
        MixedInteger(Metaheuristics.ECA(adaptive = true); options)
    else
        _check_method(model.method, multiobjective)
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

    function f_single(solution)
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

    f = if multiobjective
        mo_evaluator = MOI.Nonlinear.Evaluator(
            model.obj_model,
            MOI.Nonlinear.SparseReverseMode(),
            vars,
        )
        MOI.initialize(mo_evaluator, Symbol[])
        F = zeros(T, model.n_objectives)
        function f_multi(solution)
            x = _dict_to_values(solution, search_space, var_dist)
            MOI.eval_constraint(mo_evaluator, F, x)
            fx = objective_scale .* F # new vector: F is reused
            # multi-objective algorithms expect (f, g, h)
            has_constraints || return fx, zeros(1), zeros(1)
            gx, hx = constraints(x)
            fx, gx, hx
        end
    else
        f_single
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
