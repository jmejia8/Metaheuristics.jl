# MetaheuristicsMOIExt

[MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl) extension for
Metaheuristics. It lets you use `Metaheuristics.Optimizer` as a solver in
[JuMP](https://github.com/jump-dev/JuMP.jl).

## Installation

The extension loads automatically when `MathOptInterface` (or JuMP) is available:

```julia
julia> import Pkg; Pkg.add(["Metaheuristics", "JuMP"])
```

## Quick Start

Minimize:

$$f(x) = x_1^2 + x_2^2$$

subject to $x_1 + x_2 \geq 2$, where $x\in [-5, 5]^2$. The optimum is $f(1, 1) = 2$.

```julia
using JuMP, Metaheuristics

model = Model(Metaheuristics.Optimizer)
@variable(model, -5 <= x[1:2] <= 5)
@objective(model, Min, x[1]^2 + x[2]^2)
@constraint(model, x[1] + x[2] >= 2)

JuMP.optimize!(model)

@show objective_value(model)
@show value.(x)
```

## Mixed-Integer Constraints

Integer and binary variables are supported, as well as nonlinear constraints. The
optimum below is $y = 1$, $z = 3$, $b = 0$ with objective value $0.25$.

```julia
model = Model(Metaheuristics.Optimizer)
@variable(model, -5 <= y <= 5)
@variable(model, 0 <= z <= 10, Int)
@variable(model, b, Bin)
@objective(model, Min, (y - 1.5)^2 + (z - 3)^2 + 10b)
@constraint(model, y + z <= 4)
@constraint(model, y * z >= 1)

JuMP.optimize!(model)

@show objective_value(model)
@show value(y), value(z), value(b)
```

## Solver Options

By default, `ECA` is used within `MixedInteger`. Use the `"method"` attribute to choose
another algorithm and its options:

```julia
options = Options(seed = 1, iterations = 500)
method = MixedInteger(DE(N = 50); options)

model = Model(Metaheuristics.Optimizer)
set_attribute(model, "method", method)
set_silent(model)
```

Use `@objective(model, Max, ...)` for maximization.
