"""
Optimizer initialization and management.
"""
module OptimizerUtils

import ..Types
import ..RandomUtils
import ..Constants

"""
    initialize!(
        opt::Types.Optimizer,
        initial_pop_size::Int,
        dimension::Int,
        objective_func::Function,
        lower_bounds::Union{Float64, Vector{Float64}},
        upper_bounds::Union{Float64, Vector{Float64}};
        max_evaluations::Int = 10000 * dimension,
        func_num::Int = 0,
        func_index::Int = 0,
        seed::Union{Int, Nothing} = nothing
    )

Initialize the optimizer with problem parameters and allocate all necessary arrays.
"""
function initialize!(
    opt::Types.Optimizer,
    initial_pop_size::Int,
    dimension::Int,
    objective_func::Function,
    lower_bounds::Union{Float64, Vector{Float64}},
    upper_bounds::Union{Float64, Vector{Float64}};
    max_evaluations::Int = 10000 * dimension,
    seed::Union{Int, Nothing} = nothing
)::Nothing
    opt.dimension = dimension
    opt.objective_func = objective_func
    opt.n_evaluations = 0
    opt.max_evaluations = max_evaluations
    
    # Handle bounds: convert scalar to vector if needed
    if isa(lower_bounds, Number)
        opt.lower_bounds = fill(Float64(lower_bounds), dimension)
    else
        opt.lower_bounds = Vector{Float64}(lower_bounds)
        @assert length(opt.lower_bounds) == dimension "Lower bounds length must match dimension"
    end
    
    if isa(upper_bounds, Number)
        opt.upper_bounds = fill(Float64(upper_bounds), dimension)
    else
        opt.upper_bounds = Vector{Float64}(upper_bounds)
        @assert length(opt.upper_bounds) == dimension "Upper bounds length must match dimension"
    end
    
    # Validate bounds
    @assert all(opt.lower_bounds .< opt.upper_bounds) "Lower bounds must be less than upper bounds"
    
    # Initialize RNG
    if seed !== nothing
        opt.rng = RandomUtils.RNGState(seed)
    elseif !isdefined(opt, :rng) || opt.rng === nothing
        opt.rng = RandomUtils.RNGState()
    end
    
    # Set population parameters
    opt.n_inds_current = initial_pop_size
    opt.n_inds_front = initial_pop_size
    opt.n_inds_front_max = initial_pop_size
    opt.popul_size = initial_pop_size * 2
    opt.the_chosen_one = 0
    opt.memory_size = Constants.MEMORY_SIZE
    opt.memory_iter = 0
    opt.success_filled = 0
    opt.success_rate = 0.5
    
    # Initialize best solution
    opt.best_solution = zeros(dimension)
    opt.best_fitness = Inf
    
    # Allocate arrays
    opt.popul = zeros(opt.popul_size, dimension)
    opt.popul_front = zeros(opt.n_inds_front, dimension)
    opt.popul_temp = zeros(opt.popul_size, dimension)
    opt.fit_arr = zeros(opt.popul_size)
    opt.fit_arr_copy = zeros(opt.popul_size)
    opt.fit_arr_front = zeros(opt.n_inds_front)
    opt.trial = zeros(dimension)
    opt.temp_success_cr = zeros(opt.popul_size)
    opt.temp_success_f = zeros(opt.popul_size)
    opt.fit_delta = zeros(opt.popul_size)
    opt.memory_cr = ones(opt.memory_size)
    opt.memory_f = ones(opt.memory_size)
    opt.weights = zeros(opt.popul_size)
    opt.indices = zeros(Int, opt.popul_size)
    opt.indices2 = zeros(Int, opt.popul_size)
    opt.eb_hybrid_flag = zeros(opt.n_inds_front_max)
    opt.fit_mass = zeros(opt.n_inds_front_max)
    opt.fit_temp = zeros(opt.n_inds_front_max)
    
    # Initialize population
    for i = 1:opt.popul_size
        for j = 1:dimension
            opt.popul[i, j] = RandomUtils.random_val(opt.rng, opt.lower_bounds[j], opt.upper_bounds[j])
        end
    end
    
    opt.eb_hybrid_rate = Constants.EB_HYBRID_RATE_INIT
    return nothing
end

"""
    sort_population!(opt::Types.Optimizer, use_front::Bool = false)

Sort population by fitness. If use_front is true, sorts popul_front, otherwise sorts popul.
"""
function sort_population!(opt::Types.Optimizer, use_front::Bool = false)::Nothing
    if use_front
        n = opt.n_inds_front
        fit_arr = opt.fit_arr_front
        indices = opt.indices2
    else
        n = opt.n_inds_front
        fit_arr = opt.fit_arr
        indices = opt.indices
    end
    
    # Copy values
    for i = 1:n
        opt.fit_arr_copy[i] = fit_arr[i]
    end
    
    # Use Julia's built-in sortperm to get sorted indices
    if n > 1
        perm = sortperm(@view opt.fit_arr_copy[1:n])
        # Reorder both arrays according to permutation
        opt.fit_arr_copy[1:n] = opt.fit_arr_copy[perm]
        for i = 1:n
            indices[i] = perm[i] - 1  # Convert to 0-indexed
        end
    else
        indices[1] = 0
    end
    
    return nothing
end

"""
    initialize_front_population!(opt::Types.Optimizer)

Initialize the front population from the sorted main population.
"""
function initialize_front_population!(opt::Types.Optimizer)::Nothing
    for i = 1:opt.n_inds_front
        opt.popul_front[i, :] = opt.popul[opt.indices[i]+1, :]
        opt.fit_arr_front[i] = opt.fit_arr_copy[i]
        opt.fit_mass[i] = opt.fit_arr_front[i]
    end
    opt.pf_index = 0
    return nothing
end

"""
    truncate_population!(opt::Types.Optimizer)

Truncate population to front size if it exceeds the limit.
"""
function truncate_population!(opt::Types.Optimizer)::Nothing
    if opt.n_inds_current > opt.n_inds_front
        # Use Julia's built-in sortperm to get sorted indices
        perm = sortperm(@view opt.fit_arr[1:opt.n_inds_current])
        
        opt.n_inds_current = opt.n_inds_front
        for i = 1:opt.n_inds_current
            opt.popul_temp[i, :] = opt.popul[perm[i], :]
        end
        opt.popul[1:opt.n_inds_current, :] = opt.popul_temp[1:opt.n_inds_current, :]
        opt.fit_arr[1:opt.n_inds_current] = opt.fit_arr[perm[1:opt.n_inds_current]]
    end
    return nothing
end

end # module OptimizerUtils

