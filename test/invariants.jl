# Shared correctness checks for metaheuristic outputs (not accuracy).

using Metaheuristics
using Test

"""
    normalize_bounds(bounds)

Convert matrix bounds `[lb ub]` to `BoxConstrainedSpace` without mutating the input.
"""
function normalize_bounds(bounds)
    if bounds isa Tuple && length(bounds) == 2
        lb, ub = bounds
        eltype(lb) <: Bool && return bounds
        return BoxConstrainedSpace(lb = lb, ub = ub)
    end
    if bounds isa AbstractMatrix{Bool}
        return bounds
    end
    if bounds isa AbstractMatrix{<:Integer}
        bounds = size(bounds, 1) == 2 ? bounds : bounds'
        return BoxConstrainedSpace(lb = float.(bounds[1, :]), ub = float.(bounds[2, :]))
    end
    if bounds isa AbstractMatrix
        if size(bounds, 1) == 2
            return BoxConstrainedSpace(lb = bounds[1, :], ub = bounds[2, :])
        elseif size(bounds, 2) == 2
            return BoxConstrainedSpace(lb = bounds[:, 1], ub = bounds[:, 2])
        end
        error("Expected a 2×D or D×2 bounds matrix, got size $(size(bounds)).")
    end
    return bounds
end

function _search_space_dim(space)
    if space isa AbstractMatrix
        return size(space, 1) == 2 ? size(space, 2) : size(space, 1)
    end
    return Metaheuristics.getdim(space)
end

function _point_in_bounds(x, space)
    if space isa BoxConstrainedSpace && eltype(space.lb) <: AbstractFloat
        return all(space.lb .- 1e-6 .<= x .<= space.ub .+ 1e-6)
    end
    return x in space
end

_unwrap_solution(sol) = sol isa Metaheuristics.Bee ? sol.sol : sol

function _expected_objective(ff, sol)
    sol = _unwrap_solution(sol)
    x = get_position(sol)
    out = ff(x)
    if sol isa Metaheuristics.AbstractUnconstrainedSolution
        out isa Tuple && return out[1]
        return out
    elseif sol isa Metaheuristics.AbstractConstrainedSolution
        fx, gx, hx = out
        return fx, gx, hx
    else
        error("Unsupported solution type $(typeof(sol)) for objective_matches")
    end
end

"""
    objective_matches(ff, sol)

Check that stored objective and constraint values match re-evaluation with `ff`.
"""
function _objective_matches_constrained!(ff, sol, inner)
    fx, gx, hx = _expected_objective(ff, inner)
    @test fx ≈ fval(sol)
    @test gx ≈ gval(sol)
    @test hx ≈ hval(sol)
end

function objective_matches(ff, sol)
    inner = _unwrap_solution(sol)
    if inner isa Metaheuristics.AbstractUnconstrainedSolution
        expected = _expected_objective(ff, inner)
        @test expected ≈ fval(sol)
    elseif inner isa Metaheuristics.AbstractConstrainedSolution
        _objective_matches_constrained!(ff, sol, inner)
    else
        error("Unsupported solution type $(typeof(inner)) for objective_matches")
    end
    return true
end

function _assert_finiteness!(x, f)
    if eltype(x) <: AbstractFloat
        @test all(isfinite, x)
    end
    if f isa AbstractVector
        @test all(isfinite, f)
    else
        @test isfinite(f)
    end
end

function _assert_solution_structure!(ff, sol, space)
    x = get_position(sol)
    @test _point_in_bounds(x, space)
    @test length(x) == _search_space_dim(space)
    _assert_finiteness!(x, fval(sol))
    objective_matches(ff, sol)
end

"""
    assert_box_constrained_result!(ff, res, bounds)

Structural checks for single-objective box-constrained runs.
"""
function assert_box_constrained_result!(ff, res, bounds)
    space = normalize_bounds(bounds)
    for sol in res.population
        _assert_solution_structure!(ff, sol, space)
    end
    @test _point_in_bounds(minimizer(res), space)
    @test fval(res.best_sol) == minimum(res)
    @test minimum(res) <= minimum(fvals(res.population))
    return nothing
end

"""
    assert_multiobjective_result!(ff, res, bounds)

Structural checks for multi-objective runs.
"""
function assert_multiobjective_result!(ff, res, bounds)
    space = normalize_bounds(bounds)
    @test !isempty(res.population)
    nobj = length(fval(res.population[1]))
    for sol in res.population
        _assert_solution_structure!(ff, sol, space)
        @test length(fval(sol)) == nobj
    end
    @test _point_in_bounds(minimizer(res), space)
    pf1 = pareto_front(res)
    pf2 = pareto_front(res.population)
    @test size(pf1, 1) == size(pf2, 1) &&
          Metaheuristics.PerformanceIndicators.igd(pf1, pf2) ≈ 0.0
    return nothing
end

"""
    assert_constrained_feasibility!(res; h_tol=1e-3, g_tol=0.0)

Check constraint metadata and feasibility of solutions in `res`.
"""
function assert_constrained_feasibility!(res; h_tol=1e-3, g_tol=0.0, vio_atol=1e-2)
    for sol in res.population
        @test Metaheuristics.is_feasible(sol) == (Metaheuristics.sum_violations(sol) == 0.0)
        expected_vio = Metaheuristics.violationsSum(gval(sol), hval(sol))
        if Metaheuristics.sum_violations(sol) == 0.0
            # ε-constraint methods may mark tiny violations as feasible
            @test expected_vio < vio_atol
        else
            @test Metaheuristics.sum_violations(sol) ≈ expected_vio atol=vio_atol
        end
    end
    best = res.best_sol
    @test sum(abs.(hval(best))) < h_tol
    @test !any(gval(best) .> g_tol)
    return nothing
end

"""
    assert_constrained_result!(ff, res, bounds; h_tol, g_tol)

Bounds, stored objectives, and constraint metadata for constrained problems.
"""
function assert_constrained_result!(ff, res, bounds; h_tol=1e-3, g_tol=0.0)
    space = normalize_bounds(bounds)
    for sol in res.population
        _assert_solution_structure!(ff, sol, space)
    end
    @test _point_in_bounds(minimizer(res), space)
    @test fval(res.best_sol) == minimum(res)
    @test minimum(res) <= minimum(fvals(res.population))
    assert_constrained_feasibility!(res; h_tol=h_tol, g_tol=g_tol)
    return nothing
end

"""
    assert_minimizer_dimension!(res, dim)

Check minimizer length and finiteness without requiring convergence to bounds.
"""
function assert_minimizer_dimension!(res, dim)
    x = minimizer(res)
    @test length(x) == dim
    @test all(isfinite, x)
    @test isfinite(minimum(res))
    return nothing
end

"""
    assert_minimizer_structure!(ff, res, bounds)

Lightweight checks (minimizer only) for API / smoke tests with converged runs.
"""
function assert_minimizer_structure!(ff, res, bounds)
    space = normalize_bounds(bounds)
    sol = res.best_sol
    x = get_position(sol)
    @test _point_in_bounds(x, space)
    @test length(x) == _search_space_dim(space)
    _assert_finiteness!(x, fval(sol))
    objective_matches(ff, sol)
    @test fval(res.best_sol) == minimum(res)
    return nothing
end

"""
    assert_combinatorial_minimizer!(ff, res, bounds)

Bounds and stored objective on the reported minimizer (not the full population).
"""
function assert_combinatorial_minimizer!(ff, res, bounds)
    space = Metaheuristics.Problem(ff, bounds).search_space
    sol = res.best_sol
    x = minimizer(res)
    @test x in space
    @test length(x) == Metaheuristics.getdim(space)
    _assert_finiteness!(x, fval(sol))
    objective_matches(ff, sol)
    @test fval(res.best_sol) == minimum(res)
    return nothing
end

"""
    assert_permutation_population!(res, n)

Every population member is a valid permutation of `1:n`.
"""
function assert_permutation_population!(res, n)
    expected = collect(1:n)
    for sol in res.population
        @test sort(collect(get_position(sol))) == expected
    end
    @test sort(collect(minimizer(res))) == expected
    return nothing
end

"""
    assert_mixed_integer_result!(res, search_space)

Check MixedSpace decoding and per-variable bounds on the minimizer.
"""
function assert_mixed_integer_result!(res, search_space::MixedSpace)
    d = Metaheuristics.vec_to_dict(minimizer(res), search_space)
    @test sort(collect(keys(d))) == sort(collect(search_space.key_order))
    for k in search_space.key_order
        sub = search_space.domain[k]
        if sub isa BitArraySpace || (sub isa BoxConstrainedSpace && sub.rigid)
            @test d[k] in sub
        end
    end
    x = minimizer(res)
    @test length(x) == _search_space_dim(search_space)
    @test all(isfinite, x)
    @test all(isfinite, fval(res.best_sol))
    return nothing
end

"""
    fixed_population_size_expected(method)

Return `true` when final population size should equal `method.parameters.N`.
"""
function _base_parameters(params)
    params isa MixedInteger && return params.base
    params isa Restart && return params.alg
    return params
end

function fixed_population_size_expected(algorithm::Metaheuristics.Algorithm)
    params = algorithm.parameters
    params isa Restart && return false
    base_params = _base_parameters(params)
    base_params isa RDEx && return false
    base_params isa SHADE && return false
    if base_params isa ECA && base_params.resize_population
        return false
    end
    return hasproperty(base_params, :N) && base_params.N > 0
end

function assert_population_size!(res, algorithm::Metaheuristics.Algorithm)
    fixed_population_size_expected(algorithm) || return nothing
    @test length(res.population) == _base_parameters(algorithm.parameters).N
    return nothing
end
