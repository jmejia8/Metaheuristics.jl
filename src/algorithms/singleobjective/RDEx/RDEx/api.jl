"""
Public API for RDEx optimizer.
"""
module API

import ..Types
import ..RandomUtils
import ..OptimizerUtils
import ..OptimizeUtils

"""
    optimize(
        objective_func,
        lower_bounds,
        upper_bounds,
        max_evaluations;
        population_size=18,
        seed=nothing,
        verbose=false
    ) -> Types.OptimizationResult

Optimize an objective function using the RDEx (Random Differential Evolution with eXploration) algorithm.

# Arguments
- `objective_func::Function`: The objective function to minimize. Must accept a `Vector{Float64}` and return `Float64`.
- `lower_bounds::Union{Float64, Vector{Float64}}`: Lower bounds. If scalar, same bound for all dimensions.
- `upper_bounds::Union{Float64, Vector{Float64}}`: Upper bounds. If scalar, same bound for all dimensions.
- `max_evaluations::Int`: Maximum number of function evaluations.

# Keyword Arguments
- `population_size::Int=18`: Base population size (will be multiplied by dimension).
- `seed::Union{Int, Nothing}=nothing`: Random seed for reproducibility.
- `verbose::Bool=false`: Print progress information.

# Returns
An `OptimizationResult` with:
- `best_solution::Vector{Float64}`: Best solution found.
- `best_fitness::Float64`: Best fitness value found.
- `n_evaluations::Int`: Number of function evaluations used.
- `convergence_history::Vector{Float64}`: History of best fitness (currently empty, can be extended).

# Example
```julia
using RDEx

f(x) = sum(x.^2)  # Sphere function
result = optimize(f, -5.0, 5.0, 10000)
println("Best fitness: ", result.best_fitness)
```

# Example with vector bounds
```julia
lower = [-5.0, -10.0, -3.0]
upper = [5.0, 10.0, 3.0]
result = optimize(f, lower, upper, 10000; seed=42, verbose=true)
```
"""
function optimize(
    objective_func::Function,
    lower_bounds::Union{Float64, Vector{Float64}},
    upper_bounds::Union{Float64, Vector{Float64}},
    max_evaluations::Int;
    population_size::Int = 18,
    seed::Union{Int, Nothing} = nothing,
    verbose::Bool = false
)::Types.OptimizationResult
    # Determine dimension from bounds
    if isa(lower_bounds, Number)
        error("Cannot infer dimension from scalar bounds. Please provide vector bounds.")
    else
        dimension = length(lower_bounds)
    end
    
    # Validate bounds
    if isa(lower_bounds, Vector) && isa(upper_bounds, Vector)
        @assert length(lower_bounds) == length(upper_bounds) "Lower and upper bounds must have same length"
        @assert all(lower_bounds .< upper_bounds) "Lower bounds must be less than upper bounds"
    end
    
    # Calculate initial population size
    initial_pop_size = population_size * dimension
    
    # Create optimizer with minimal initialization
    # The initialize! function will properly allocate all arrays
    opt = Types.Optimizer(
        objective_func,           # objective_func
        Float64[],                # lower_bounds (will be set in initialize!)
        Float64[],                # upper_bounds (will be set in initialize!)
        0,                        # dimension (will be set in initialize!)
        0,                        # n_evaluations
        max_evaluations,          # max_evaluations
        RandomUtils.RNGState(),    # rng (will be set in initialize! if seed provided)
        5, 0, 0, 0, 0,           # memory parameters
        0, 0, 0, 0, 0,            # population dimensions
        0, 0,                     # selection indices
        Inf,                      # best_fitness
        Float64[],                # best_solution
        0.0, 0.0,                 # adaptation parameters (success_rate, eb_hybrid_rate)
        # Arrays - initialize! will replace these
        zeros(1,1), zeros(1,1), zeros(1,1),  # popul, popul_front, popul_temp
        zeros(1), zeros(1), zeros(1), zeros(1),  # fit_arr, fit_arr_copy, fit_arr_front, trial
        zeros(1), zeros(1), zeros(1), zeros(1), zeros(1),  # temp_success_cr, temp_success_f, memory_cr, memory_f, fit_delta
        zeros(1),  # weights
        zeros(Int,1), zeros(Int,1),  # indices, indices2
        zeros(1), zeros(1), zeros(1)  # fit_mass, fit_temp, eb_hybrid_flag
    )
    
    # Initialize optimizer (allocates arrays, sets bounds, sets up population)
    OptimizerUtils.initialize!(
        opt, initial_pop_size, dimension, objective_func, 
        lower_bounds, upper_bounds;
        max_evaluations=max_evaluations,
        seed=seed
    )
    
    if verbose
        println("Starting RDEx optimization...")
        println("  Dimension: $dimension")
        println("  Population size: $(opt.popul_size)")
        println("  Max evaluations: $max_evaluations")
    end
    
    # Run optimization
    OptimizeUtils.optimize_cycle!(opt)
    
    if verbose
        println("Optimization completed.")
        println("  Best fitness: $(opt.best_fitness)")
        println("  Evaluations used: $(opt.n_evaluations)")
    end
    
    # Return results
    return Types.OptimizationResult(
        copy(opt.best_solution),
        opt.best_fitness,
        opt.n_evaluations,
        Float64[]  # convergence_history - could be populated if tracking is added
    )
end

end # module API

