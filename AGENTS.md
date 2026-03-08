# AGENTS.md - Metaheuristics.jl Development Guide

This file provides guidelines for agentic coding agents working on the Metaheuristics.jl repository.

## Project Overview

Metaheuristics.jl is a Julia package for optimization algorithms (DE, PSO, ABC, WOA, ECA, NSGA2, etc.). It supports single-objective, multi-objective, and combinatorial optimization.

---

## Build, Lint, and Test Commands

### Running Tests

```bash
# Run all tests
julia --project -e 'include("test/runtests.jl")'

# Run a specific test file
julia --project -e 'include("test/box-constrained.jl")'

# Run a specific test set from REPL
julia --project
julia> using Test
julia> @testset "Box-constrained" begin include("test/box-constrained.jl") end
```

### Test Structure

- Tests are in `test/` directory
- Main entry: `test/runtests.jl` (uses Aqua.jl for code quality checks)
- Individual test files: `box-constrained.jl`, `constrained.jl`, `multi-objective.jl`, etc.
- Test helper function pattern: see `test/box-constrained.jl:9` for `run_methods()`

### Building/Precompiling

```bash
# Build the package
julia --project -e 'using Pkg; Pkg.build("Metaheuristics")'

# Precompile
julia --project -e 'using Metaheuristics; @time using Metaheuristics'
```

---

## Code Style Guidelines

### Naming Conventions

- **Types**: CamelCase (e.g., `DE`, `Options`, `Information`, `State`)
- **Abstract types**: Prefix with `Abstract` (e.g., `AbstractSolution`, `AbstractAlgorithm`)
- **Functions**: lowercase with underscores where appropriate (e.g., `update_state!`, `is_better`)
- **Constants**: UPPERCASE (e.g., `TerminationStatusCode`)
- **Modules/Algorithms**: Capitalized (e.g., `NSGA2`, `PSO`)

### File Organization

```
src/
├── Metaheuristics.jl          # Main module, exports
├── core/                      # Abstract types, core structures
├── solutions/                 # Solution/Individual types
├── operators/                 # Evolutionary operators
├── termination/               # Termination criteria
├── common/                    # Shared utilities
├── algorithms/
│   ├── singleobjective/       # DE, PSO, ABC, WOA, etc.
│   ├── multiobjective/        # NSGA2, SMS_EMOA, etc.
│   └── combinatorial/         # LocalSearch, VNS, etc.
└── ...
```

### Imports

- Use `using` for imports that are frequently used
- Use qualified imports for ambiguous names: `using Random: rand`
- Group stdlib imports first, then external packages

```julia
using Test
using Random
using LinearAlgebra
using Metaheuristics
```

### Documentation

- Use triple-quoted strings for docstrings
- Include examples in docstrings using `julia> ...` blocks
- Document public API functions thoroughly

```julia
"""
    DE(; N=0, F=1.0, CR=0.5, strategy=:rand1, ...)

Parameters for Differential Evolution algorithm.

# Example
```jldoctest
julia> f(x) = sum(x.^2)
julia> optimize(f, [-1 -1 -1; 1 1 1.0], DE())
```
"""
function DE(; ...) end
```

### Type Annotations

- Provide type annotations for struct fields
- Use abstract types for interfaces (e.g., `AbstractProblem`)
- Use `Union{T, Nothing}` for optional fields
- Be consistent with Float64 vs Int types

```julia
mutable struct DE <: AbstractDifferentialEvolution
    N::Int
    F::Float64
    CR::Float64
    strategy::Symbol
end
```

### Function Design

- Use keyword arguments for algorithm parameters
- Use mutating functions with `!` suffix where appropriate
- Follow Julia's multiple dispatch pattern

```julia
# Constructor pattern
function DE(; N::Int=0, F=0.7, CR=0.5, kargs...)
    parameters = DE(N, F, CR)
    Algorithm(parameters; kargs...)
end

# Update function pattern
function update_state!(status, parameters::AbstractDifferentialEvolution, 
                      problem::AbstractProblem, information::Information, 
                      options::Options, args...; kargs...)
    # implementation
end
```

### Error Handling

- Use `@assert` for debugging/invariant checks (not for user input validation)
- Use `@warn` for non-fatal warnings with `options.debug` check
- Use `throw()` for fatal errors

```julia
@assert !isempty(population)
options.debug && @warn("CR should be from interval [0,1]; set to default value 0.5")
```

### Algorithm Interface

New algorithms must implement:

1. **Constructor**: `Algorithm(params; information=Information(), options=Options())`
2. **`initialize!`**: Initialize population and state
3. **`update_state!`**: Main iteration logic
4. **`final_stage!`**: Cleanup, timing
5. **`is_better`**: Solution comparison

Example structure:

```julia
# Parameters struct
mutable struct MyAlgo <: AbstractParameters
    N::Int
    param::Float64
end

# Constructor
function MyAlgo(; N=0, param=0.5, kargs...)
    Algorithm(MyAlgo(N, param); kargs...)
end

# Required methods
function initialize!(status, params::MyAlgo, problem, information, options, args...; kargs...)
    # Initialize population
end

function update_state!(status, params::MyAlgo, problem, information, options, args...; kargs...)
    # Main loop
end

function final_stage!(status, params::MyAlgo, problem, information, options, args...; kargs...)
    status.final_time = time()
end

is_better(a, b, params::MyAlgo) = is_better(a, b)
```

### Testing New Algorithms

- Add tests in appropriate file in `test/`
- Test basic functionality, convergence, edge cases
- Use `Metaheuristics.TestProblems.get_problem()` for benchmark functions

```julia
@testset "MyAlgo" begin
    f(x) = sum(x.^2)
    bounds = [-1 -1; 1 1]
    result = optimize(f, bounds, MyAlgo())
    @test minimum(result) < 1e-4
end
```

### Performance Considerations

- Pre-allocate vectors when possible
- Use views (`@views`) to avoid copying
- Minimize allocations in hot loops
- Use `rand(options.rng)` for reproducible random numbers

---

## Additional Notes

- The package uses SnoopPrecompile for faster load times
- Aqua.jl is used for code quality checks (ambiguities, etc.)
- Document any new public API exports in `src/Metaheuristics.jl`
- Follow the existing code patterns when extending functionality
