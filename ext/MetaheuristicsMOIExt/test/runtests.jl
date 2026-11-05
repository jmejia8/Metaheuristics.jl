# Tests for the MathOptInterface extension (MetaheuristicsMOIExt).
#
# Kept apart from test/ so `Pkg.test()` is unaffected. Run from the repo root:
#
#   julia --project=ext/MetaheuristicsMOIExt/test -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
#   julia --project=ext/MetaheuristicsMOIExt/test ext/MetaheuristicsMOIExt/test/runtests.jl

using Test
using Metaheuristics
import MathOptInterface as MOI

const Ext = Base.get_extension(Metaheuristics, :MetaheuristicsMOIExt)

# Small, seeded solver so tests are fast and reproducible.
function _method(; seed = 1, iterations = 300)
    options = Options(; seed, iterations, f_calls_limit = Inf)
    MixedInteger(ECA(N = 30); options)
end

function _new_optimizer(; kwargs...)
    opt = Ext.Optimizer()
    MOI.set(opt, MOI.Silent(), true)
    MOI.set(opt, MOI.RawOptimizerAttribute("method"), _method(; kwargs...))
    opt
end

_sum_of_squares(xs, c = 0.0) = MOI.ScalarNonlinearFunction(
    :+,
    Any[MOI.ScalarNonlinearFunction(:^, Any[x - c, 2]) for x in xs]
)

@testset "MetaheuristicsMOIExt" begin

    @testset "Extension loading" begin
        @test Ext !== nothing
        @test isdefined(Metaheuristics, :Optimizer)
        @test Metaheuristics.Optimizer === Ext.Optimizer
    end

    @testset "Attributes" begin
        opt = Ext.Optimizer()
        @test opt isa MOI.AbstractOptimizer
        @test opt isa Ext.Optimizer{Float64}
        @test MOI.is_empty(opt)
        @test MOI.get(opt, MOI.SolverName()) == "Metaheuristics"
        @test MOI.get(opt, MOI.TerminationStatus()) == MOI.OPTIMIZE_NOT_CALLED
        @test MOI.get(opt, MOI.ResultCount()) == 0
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.NO_SOLUTION
        @test MOI.get(opt, MOI.DualStatus()) == MOI.NO_SOLUTION

        @test MOI.supports(opt, MOI.Silent())
        MOI.set(opt, MOI.Silent(), true)
        @test MOI.get(opt, MOI.Silent())

        @test MOI.supports(opt, MOI.TimeLimitSec())
        @test MOI.get(opt, MOI.TimeLimitSec()) === nothing
        MOI.set(opt, MOI.TimeLimitSec(), 5)
        @test MOI.get(opt, MOI.TimeLimitSec()) == 5.0
        MOI.set(opt, MOI.TimeLimitSec(), nothing)
        @test MOI.get(opt, MOI.TimeLimitSec()) === nothing

        @test MOI.supports(opt, MOI.RawOptimizerAttribute("method"))
        @test MOI.supports(opt, MOI.RawOptimizerAttribute("iterations"))
        @test !MOI.supports(opt, MOI.RawOptimizerAttribute("not_an_option"))
        MOI.set(opt, MOI.RawOptimizerAttribute("iterations"), 10)
        @test MOI.get(opt, MOI.RawOptimizerAttribute("iterations")) == 10
        @test_throws ErrorException MOI.get(opt, MOI.RawOptimizerAttribute("seed"))

        # "method" getter: nothing until set, then the exact object that was set
        @test MOI.get(opt, MOI.RawOptimizerAttribute("method")) === nothing
        method = _method()
        MOI.set(opt, MOI.RawOptimizerAttribute("method"), method)
        @test MOI.get(opt, MOI.RawOptimizerAttribute("method")) === method
    end

    @testset "Supported constraints" begin
        opt = Ext.Optimizer()
        for S in (
                MOI.LessThan{Float64},
                MOI.GreaterThan{Float64},
                MOI.EqualTo{Float64},
                MOI.Interval{Float64},
                MOI.ZeroOne,
                MOI.Integer,
            )
            @test MOI.supports_constraint(opt, MOI.VariableIndex, S)
        end
        for S in (MOI.LessThan{Float64}, MOI.GreaterThan{Float64}, MOI.EqualTo{Float64})
            @test MOI.supports_constraint(opt, MOI.ScalarNonlinearFunction, S)
        end
        @test !MOI.supports(opt, MOI.VariablePrimalStart(), MOI.VariableIndex)
    end

    @testset "Model building and empty!" begin
        opt = _new_optimizer()
        x = MOI.add_variables(opt, 2)
        @test MOI.is_valid(opt, x[1])
        MOI.add_constraint(opt, x[1], MOI.Interval(-1.0, 1.0))
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x),
        )
        @test !MOI.is_empty(opt)
        MOI.empty!(opt)
        @test MOI.is_empty(opt)
    end

    @testset "optimize! without objective is a feasibility problem" begin
        opt = _new_optimizer()
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-1.0, 1.0))
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test -1 <= MOI.get(opt, MOI.VariablePrimal(), x) <= 1
        @test MOI.get(opt, MOI.ObjectiveValue()) == 0
    end

    @testset "Continuous, box constrained (min)" begin
        opt = _new_optimizer()
        x = MOI.add_variables(opt, 3)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
        end
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x, 1.0),
        )
        MOI.optimize!(opt)

        @test MOI.get(opt, MOI.ResultCount()) == 1
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(opt, MOI.TerminationStatus()) in
            (MOI.LOCALLY_SOLVED, MOI.ITERATION_LIMIT, MOI.OTHER_LIMIT)
        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ 0.0 atol = 1.0e-3
        @test MOI.get.(opt, MOI.VariablePrimal(), x) ≈ ones(3) atol = 1.0e-2
        @test MOI.get(opt, MOI.SolveTimeSec()) >= 0
    end

    @testset "Maximization sense" begin
        opt = _new_optimizer()
        x = MOI.add_variables(opt, 2)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-2.0, 2.0))
        end
        # max 4 - x1^2 - x2^2  ->  optimum 4 at (0, 0)
        obj = MOI.ScalarNonlinearFunction(:-, Any[4.0, _sum_of_squares(x)])
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(), obj)
        MOI.optimize!(opt)

        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ 4.0 atol = 1.0e-3
        @test MOI.get.(opt, MOI.VariablePrimal(), x) ≈ zeros(2) atol = 1.0e-2
    end

    @testset "Unbounded variables (non-rigid space)" begin
        opt = _new_optimizer(iterations = 500)
        x = MOI.add_variables(opt, 2)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x, 3.0),
        )
        MOI.optimize!(opt)

        @test MOI.get.(opt, MOI.VariablePrimal(), x) ≈ fill(3.0, 2) atol = 1.0e-2
    end

    @testset "Nonlinear constraint" begin
        opt = _new_optimizer(iterations = 500)
        x = MOI.add_variables(opt, 2)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
        end
        # min x1^2 + x2^2  s.t.  x1 + x2 >= 2  ->  optimum 2 at (1, 1)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x),
        )
        g = MOI.ScalarNonlinearFunction(:+, Any[x[1], x[2]])
        ci = MOI.add_constraint(opt, g, MOI.GreaterThan(2.0))
        @test ci isa
            MOI.ConstraintIndex{MOI.ScalarNonlinearFunction, MOI.GreaterThan{Float64}}
        MOI.optimize!(opt)

        xv = MOI.get.(opt, MOI.VariablePrimal(), x)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test sum(xv) >= 2 - 1.0e-6
        # loose: ECA's convergence test may stop early on constrained problems
        # (core library; this seed gives 2.05, other seeds stop farther away)
        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ 2.0 atol = 1.0e-1
    end

    @testset "Mixed integer and binary variables" begin
        opt = _new_optimizer(iterations = 500)
        y = MOI.add_variable(opt) # continuous
        z = MOI.add_variable(opt) # integer
        b = MOI.add_variable(opt) # binary
        MOI.add_constraint(opt, y, MOI.Interval(-5.0, 5.0))
        MOI.add_constraint(opt, z, MOI.Interval(-5.0, 5.0))
        MOI.add_constraint(opt, z, MOI.Integer())
        MOI.add_constraint(opt, b, MOI.ZeroOne())
        # min (y - 0.5)^2 + (z - 2)^2 + (b - 1)^2  ->  (0.5, 2, 1)
        obj = MOI.ScalarNonlinearFunction(
            :+,
            Any[
                MOI.ScalarNonlinearFunction(:^, Any[y - 0.5, 2]),
                MOI.ScalarNonlinearFunction(:^, Any[z - 2.0, 2]),
                MOI.ScalarNonlinearFunction(:^, Any[b - 1.0, 2]),
            ],
        )
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(), obj)
        MOI.optimize!(opt)

        yv, zv, bv = MOI.get.(opt, MOI.VariablePrimal(), (y, z, b))
        @test yv ≈ 0.5 atol = 1.0e-2
        @test isinteger(zv) && zv == 2
        @test bv in (0, 1) && bv == 1
    end

    # Regression tests for the review of 2026-11-05 (items 1 to 3).

    @testset "Re-solve after objective change uses a fresh population" begin
        opt = _new_optimizer(iterations = 200)
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-5.0, 5.0))
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        F = MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}()
        MOI.set(opt, F, _sum_of_squares([x], 1.0))
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.VariablePrimal(), x) ≈ 1.0 atol = 1.0e-2

        MOI.set(opt, F, _sum_of_squares([x], -3.0))
        MOI.optimize!(opt)
        xv = MOI.get(opt, MOI.VariablePrimal(), x)
        @test xv ≈ -3.0 atol = 1.0e-2
        # reported objective must match the current objective at the minimizer
        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ (xv + 3)^2 atol = 1.0e-8
    end

    @testset "Re-solve after adding a variable" begin
        opt = _new_optimizer(iterations = 200)
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-5.0, 5.0))
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        F = MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}()
        MOI.set(opt, F, _sum_of_squares([x]))
        MOI.optimize!(opt)

        y = MOI.add_variable(opt)
        MOI.add_constraint(opt, y, MOI.Interval(-5.0, 5.0))
        MOI.set(opt, F, _sum_of_squares([x, y], 1.0))
        MOI.optimize!(opt)
        @test MOI.get.(opt, MOI.VariablePrimal(), [x, y]) ≈ ones(2) atol = 1.0e-2
    end

    @testset "Method object is not mutated and can be shared" begin
        method = _method(iterations = 100)
        sols = map((2, 3)) do n
            opt = Ext.Optimizer()
            MOI.set(opt, MOI.Silent(), true)
            MOI.set(opt, MOI.RawOptimizerAttribute("method"), method)
            x = MOI.add_variables(opt, n)
            for xi in x
                MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
            end
            MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(
                opt,
                MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
                _sum_of_squares(x),
            )
            MOI.optimize!(opt)
            @test MOI.get(opt, MOI.RawOptimizerAttribute("method")) === method
            MOI.get.(opt, MOI.VariablePrimal(), x)
        end
        @test length(sols[2]) == 3
        @test sols[2] ≈ zeros(3) atol = 1.0e-2
        @test isnothing(method.status.best_sol)
        @test isempty(method.status.population)
    end

    @testset "Re-solving the same model is reproducible" begin
        opt = _new_optimizer(iterations = 50, seed = 7)
        x = MOI.add_variables(opt, 3)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
        end
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x, 1.0),
        )
        MOI.optimize!(opt)
        f1 = MOI.get(opt, MOI.ObjectiveValue())
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.ObjectiveValue()) == f1
    end

    @testset "Semi-bounded variables without nonlinear constraints" begin
        # min x1 - x2  s.t.  x1 >= 1, x2 <= -2  ->  (1, -2), unbounded without the bounds
        for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
            opt = _new_optimizer(iterations = 500)
            x = MOI.add_variables(opt, 2)
            MOI.add_constraint(opt, x[1], MOI.GreaterThan(1.0))
            MOI.add_constraint(opt, x[2], MOI.LessThan(-2.0))
            s = sense == MOI.MIN_SENSE ? 1.0 : -1.0
            obj = MOI.ScalarAffineFunction(
                [MOI.ScalarAffineTerm(s, x[1]), MOI.ScalarAffineTerm(-s, x[2])],
                0.0,
            )
            MOI.set(opt, MOI.ObjectiveSense(), sense)
            MOI.set(opt, MOI.ObjectiveFunction{typeof(obj)}(), obj)
            MOI.optimize!(opt)
            @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
            xv = MOI.get.(opt, MOI.VariablePrimal(), x)
            @test xv[1] >= 1.0 && xv[2] <= -2.0
            @test xv ≈ [1.0, -2.0] atol = 1.0e-2
        end
    end

    @testset "Bounds on binary variables are enforced" begin
        # fixed binary: b = 1, min b  ->  b = 1
        opt = _new_optimizer()
        b = MOI.add_variable(opt)
        MOI.add_constraint(opt, b, MOI.EqualTo(1.0))
        MOI.add_constraint(opt, b, MOI.ZeroOne())
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.VariableIndex}(), b)
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.VariablePrimal(), b) == 1
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT

        # b <= 0, max b  ->  b = 0
        opt = _new_optimizer()
        b = MOI.add_variable(opt)
        MOI.add_constraint(opt, b, MOI.LessThan(0.0))
        MOI.add_constraint(opt, b, MOI.ZeroOne())
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MAX_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.VariableIndex}(), b)
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.VariablePrimal(), b) == 0
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
    end

    @testset "ConstraintPrimal of variable constraints, mixed types" begin
        opt = _new_optimizer()
        y = MOI.add_variable(opt) # continuous, MOI index 1
        b = MOI.add_variable(opt) # binary, MOI index 2 (first in MixedSpace order)
        z = MOI.add_variable(opt) # integer, MOI index 3
        c_ub = MOI.add_constraint(opt, y, MOI.LessThan(5.0))
        c_lb = MOI.add_constraint(opt, y, MOI.GreaterThan(0.0))
        c_bin = MOI.add_constraint(opt, b, MOI.ZeroOne())
        c_int_box = MOI.add_constraint(opt, z, MOI.Interval(-3.0, 3.0))
        c_int = MOI.add_constraint(opt, z, MOI.Integer())
        # min (y - 4)^2 + (1 - b) + (z + 2)^2  ->  (4, 1, -2)
        obj = MOI.ScalarNonlinearFunction(
            :+,
            Any[
                MOI.ScalarNonlinearFunction(:^, Any[y - 4.0, 2]),
                MOI.ScalarNonlinearFunction(:-, Any[1.0, b]),
                MOI.ScalarNonlinearFunction(:^, Any[z + 2.0, 2]),
            ],
        )
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(), obj)
        MOI.optimize!(opt)

        yv, bv, zv = MOI.get.(opt, MOI.VariablePrimal(), (y, b, z))
        @test yv ≈ 4.0 atol = 1.0e-2
        @test bv == 1
        @test zv == -2
        @test MOI.get(opt, MOI.ConstraintPrimal(), c_ub) == yv
        @test MOI.get(opt, MOI.ConstraintPrimal(), c_lb) == yv
        @test MOI.get(opt, MOI.ConstraintPrimal(), c_bin) == bv
        @test MOI.get(opt, MOI.ConstraintPrimal(), c_int_box) == zv
        @test MOI.get(opt, MOI.ConstraintPrimal(), c_int) == zv
    end

    # Regression tests for the review of 2026-11-05 (items 4 to 10).

    # min sum(x.^2) on [-5, 5]^n, used by the option tests
    function _box_model(opt, n)
        x = MOI.add_variables(opt, n)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
        end
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x),
        )
        x
    end

    @testset "Raw attributes and TimeLimitSec are applied" begin
        # "iterations" overrides the method options, which stay untouched
        opt = _new_optimizer(iterations = 300)
        method = MOI.get(opt, MOI.RawOptimizerAttribute("method"))
        _box_model(opt, 3)
        MOI.set(opt, MOI.RawOptimizerAttribute("iterations"), 5)
        MOI.optimize!(opt)
        @test opt.results.iteration == 5
        @test MOI.get(opt, MOI.TerminationStatus()) == MOI.ITERATION_LIMIT
        @test method.options.iterations == 300

        # time limit, also with the default method
        opt = Ext.Optimizer()
        MOI.set(opt, MOI.Silent(), true)
        x = MOI.add_variables(opt, 30)
        for xi in x
            MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
        end
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        obj = MOI.ScalarNonlinearFunction(
            :+,
            Any[MOI.ScalarNonlinearFunction(:^, Any[xi, 4]) for xi in x]
        )
        MOI.set(opt, MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(), obj)
        MOI.set(opt, MOI.TimeLimitSec(), 0.2)
        MOI.set(opt, MOI.RawOptimizerAttribute("f_tol"), 0.0)
        MOI.set(opt, MOI.RawOptimizerAttribute("x_tol"), 0.0)
        MOI.set(opt, MOI.RawOptimizerAttribute("f_tol_rel"), 0.0)
        t = @elapsed MOI.optimize!(opt)
        @test MOI.get(opt, MOI.TerminationStatus()) == MOI.TIME_LIMIT
        @test MOI.get(opt, MOI.SolveTimeSec()) < 1.0
        @test t < 2.0

        # "seed" gives reproducible results with the default method
        objs = map(1:2) do _
            opt = Ext.Optimizer()
            MOI.set(opt, MOI.Silent(), true)
            MOI.set(opt, MOI.RawOptimizerAttribute("seed"), 3)
            MOI.set(opt, MOI.RawOptimizerAttribute("iterations"), 20)
            _box_model(opt, 3)
            MOI.optimize!(opt)
            MOI.get(opt, MOI.ObjectiveValue())
        end
        @test objs[1] == objs[2]

        # unsupported names and wrong value types fail when set
        opt = Ext.Optimizer()
        @test_throws MOI.UnsupportedAttribute MOI.set(
            opt,
            MOI.RawOptimizerAttribute("not_an_option"),
            1,
        )
        @test_throws InexactError MOI.set(opt, MOI.RawOptimizerAttribute("iterations"), 2.5)
        @test MOI.get(opt, MOI.SolveTimeSec()) == 0.0
    end

    @testset "Infeasible result is not LOCALLY_SOLVED" begin
        # x in [-1, 1], x^2 >= 4: no feasible point
        opt = _new_optimizer()
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-1.0, 1.0))
        MOI.add_constraint(
            opt,
            MOI.ScalarNonlinearFunction(:^, Any[x, 2]),
            MOI.GreaterThan(4.0),
        )
        MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
        MOI.set(opt, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.INFEASIBLE_POINT
        @test MOI.get(opt, MOI.TerminationStatus()) in
            (MOI.LOCALLY_INFEASIBLE, MOI.ITERATION_LIMIT, MOI.OTHER_LIMIT)
    end

    @testset "Integrality does not depend on constraint order" begin
        for integer_first in (true, false)
            opt = _new_optimizer()
            z = MOI.add_variable(opt)
            integer_first && MOI.add_constraint(opt, z, MOI.Integer())
            MOI.add_constraint(opt, z, MOI.Interval(0.0, 10.0))
            integer_first || MOI.add_constraint(opt, z, MOI.Integer())
            MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(
                opt,
                MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
                MOI.ScalarNonlinearFunction(:^, Any[1.0 * z - 2.4, 2]),
            )
            MOI.optimize!(opt)
            @test MOI.get(opt, MOI.VariablePrimal(), z) == 2
        end
    end

    @testset "Integer variables with fractional or huge bounds" begin
        _int_model(lb, ub, sense) = begin
            opt = _new_optimizer()
            z = MOI.add_variable(opt)
            MOI.add_constraint(opt, z, MOI.Interval(lb, ub))
            MOI.add_constraint(opt, z, MOI.Integer())
            MOI.set(opt, MOI.ObjectiveSense(), sense)
            MOI.set(opt, MOI.ObjectiveFunction{MOI.VariableIndex}(), z)
            MOI.optimize!(opt)
            opt, z
        end

        opt, z = _int_model(0.5, 3.5, MOI.MIN_SENSE)
        @test MOI.get(opt, MOI.VariablePrimal(), z) == 1
        opt, z = _int_model(0.5, 3.5, MOI.MAX_SENSE)
        @test MOI.get(opt, MOI.VariablePrimal(), z) == 3

        # huge finite bounds do not fit in Int
        opt, z = _int_model(-1.0e20, 2.5, MOI.MAX_SENSE)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test MOI.get(opt, MOI.VariablePrimal(), z) == 2

        # no integer in [0.2, 0.8]: infeasible, without calling the solver
        opt, z = _int_model(0.2, 0.8, MOI.MIN_SENSE)
        @test MOI.get(opt, MOI.TerminationStatus()) == MOI.INFEASIBLE
        @test MOI.get(opt, MOI.ResultCount()) == 0
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.NO_SOLUTION
        @test occursin("integer variable", MOI.get(opt, MOI.RawStatusString()))
    end

    @testset "Equality constraints use h_tol" begin
        for h_tol in (nothing, 1.0e-3)
            opt = _new_optimizer(iterations = 500)
            isnothing(h_tol) || MOI.set(opt, MOI.RawOptimizerAttribute("h_tol"), h_tol)
            x = MOI.add_variables(opt, 2)
            for xi in x
                MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
            end
            # min x1^2 + x2^2  s.t.  x1 + 2x2 == 2  ->  (0.4, 0.8)
            MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            MOI.set(
                opt,
                MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
                _sum_of_squares(x),
            )
            h = MOI.ScalarNonlinearFunction(:+, Any[x[1], 2.0 * x[2]])
            MOI.add_constraint(opt, h, MOI.EqualTo(2.0))
            MOI.optimize!(opt)
            xv = MOI.get.(opt, MOI.VariablePrimal(), x)
            tol = something(h_tol, Ext.DEFAULT_H_TOL)
            @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
            @test abs(xv[1] + 2xv[2] - 2) <= tol
        end
    end

    @testset "Feasibility problems" begin
        # constraints and no objective: x in [-1, 1], x^2 >= 0.25
        opt = _new_optimizer()
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-1.0, 1.0))
        MOI.add_constraint(
            opt,
            MOI.ScalarNonlinearFunction(:^, Any[x, 2]),
            MOI.GreaterThan(0.25),
        )
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test abs(MOI.get(opt, MOI.VariablePrimal(), x)) >= 0.5

        # FEASIBILITY_SENSE ignores an objective that was set
        opt = _new_optimizer()
        x = MOI.add_variable(opt)
        MOI.add_constraint(opt, x, MOI.Interval(-1.0, 1.0))
        MOI.set(opt, MOI.ObjectiveFunction{MOI.VariableIndex}(), x)
        MOI.set(opt, MOI.ObjectiveSense(), MOI.FEASIBILITY_SENSE)
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.ObjectiveValue()) == 0

        # nothing at all (no objective, no constraints): no error
        opt = Ext.Optimizer()
        MOI.set(opt, MOI.Silent(), true)
        MOI.set(opt, MOI.RawOptimizerAttribute("iterations"), 10)
        MOI.add_variable(opt)
        MOI.optimize!(opt)
        @test MOI.get(opt, MOI.ResultCount()) == 1
    end

    @testset "Legacy NLPBlock" begin
        _block(nlp, x) = MOI.NLPBlockData(MOI.Nonlinear.Evaluator(
            nlp,
            MOI.Nonlinear.SparseReverseMode(),
            x,
        ))
        _new(n) = begin
            opt = _new_optimizer(iterations = 500)
            x = MOI.add_variables(opt, n)
            for xi in x
                MOI.add_constraint(opt, xi, MOI.Interval(-5.0, 5.0))
            end
            MOI.set(opt, MOI.ObjectiveSense(), MOI.MIN_SENSE)
            opt, x
        end

        @test MOI.supports(Ext.Optimizer(), MOI.NLPBlock())

        # objective and constraint in the block (@NLobjective, @NLconstraint):
        # min (x1 - 3)^2 + (x2 + 1)^2  s.t.  x1 + x2 <= 1  ->  (2.5, -1.5)
        opt, x = _new(2)
        nlp = MOI.Nonlinear.Model()
        MOI.Nonlinear.set_objective(nlp, :(($(x[1]) - 3)^2 + ($(x[2]) + 1)^2))
        MOI.Nonlinear.add_constraint(nlp, :($(x[1]) + $(x[2])), MOI.LessThan(1.0))
        block = _block(nlp, x)
        MOI.set(opt, MOI.NLPBlock(), block)
        @test MOI.get(opt, MOI.NLPBlock()) === block
        MOI.optimize!(opt)
        xv = MOI.get.(opt, MOI.VariablePrimal(), x)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test sum(xv) <= 1 + 1.0e-8
        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ 0.5 atol = 1.0e-1

        # constraint in the block, objective from the new interface
        # (@NLconstraint with @objective): min x1^2 + x2^2  s.t.  x1 + x2 >= 2
        opt, x = _new(2)
        nlp = MOI.Nonlinear.Model()
        MOI.Nonlinear.add_constraint(nlp, :($(x[1]) + $(x[2])), MOI.GreaterThan(2.0))
        MOI.set(opt, MOI.NLPBlock(), _block(nlp, x))
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x),
        )
        MOI.optimize!(opt)
        xv = MOI.get.(opt, MOI.VariablePrimal(), x)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test sum(xv) >= 2 - 1.0e-8
        @test MOI.get(opt, MOI.ObjectiveValue()) ≈ 2.0 atol = 1.0e-1

        # block constraints are combined with new-interface constraints
        opt, x = _new(2)
        nlp = MOI.Nonlinear.Model()
        MOI.Nonlinear.add_constraint(nlp, :($(x[1])), MOI.GreaterThan(1.0))
        MOI.set(opt, MOI.NLPBlock(), _block(nlp, x))
        MOI.add_constraint(
            opt,
            MOI.ScalarNonlinearFunction(:+, Any[x[2]]),
            MOI.GreaterThan(2.0),
        )
        MOI.set(
            opt,
            MOI.ObjectiveFunction{MOI.ScalarNonlinearFunction}(),
            _sum_of_squares(x),
        )
        MOI.optimize!(opt)
        xv = MOI.get.(opt, MOI.VariablePrimal(), x)
        @test MOI.get(opt, MOI.PrimalStatus()) == MOI.FEASIBLE_POINT
        @test xv[1] >= 1 && xv[2] >= 2

        MOI.empty!(opt)
        @test isnothing(MOI.get(opt, MOI.NLPBlock()))
        @test MOI.is_empty(opt)
    end
end
