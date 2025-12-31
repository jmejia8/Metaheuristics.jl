"""
Metaheuristics.jl API integration for RDEx.

This module provides the interface to use RDEx as a standard Metaheuristics algorithm.
"""
module MetaheuristicsAPI

using Metaheuristics
using Random
using Distributions
import Metaheuristics: initialize!, update_state!, final_stage!
import Metaheuristics: AbstractParameters, gen_initial_state, Algorithm, create_solution, create_solutions, is_better
import Metaheuristics: get_position, positions, fvals, get_best, getdim

import ..SelectionUtils
import ..MutationUtils
import ..AdaptationUtils
import ..Constants

"""
    RDEx

Parameters for the RDEx algorithm compatible with Metaheuristics.jl API.

# Fields
- `N::Int`: Population size (calculated as population_size * dimension in initialize!)
- `population_size::Int`: Base population size multiplier (default 18)
- `memory_size::Int`: Memory size for parameter adaptation (default 5)
- `n_inds_front::Int`: Current front population size
- `n_inds_front_max::Int`: Maximum front population size
- `success_rate::Float64`: Success rate for adaptation
- `eb_hybrid_rate::Float64`: Elite-based hybrid rate
- `memory_iter::Int`: Current memory iteration index
- `success_filled::Int`: Number of successful trials in current generation
- `memory_cr::Vector{Float64}`: Memory for crossover rates
- `memory_f::Vector{Float64}`: Memory for scaling factors
- `temp_success_cr::Vector{Float64}`: Temporary storage for successful CR values
- `temp_success_f::Vector{Float64}`: Temporary storage for successful F values
- `fit_delta::Vector{Float64}`: Fitness differences for weighted mean calculation
- `weights::Vector{Float64}`: Weights for adaptation
- `fit_mass::Vector{Float64}`: Fitness mass for hybrid parameter update
- `fit_temp::Vector{Float64}`: Temporary fitness storage
- `eb_hybrid_flag::Vector{Float64}`: Flags indicating elite-based mutation usage
- `trial::Vector{Float64}`: Trial vector for mutation
- `indices::Vector{Int}`: Indices for sorting
- `indices2::Vector{Int}`: Secondary indices for sorting
"""
mutable struct RDEx <: AbstractParameters
    N::Int  # Population size (set in initialize!)
    population_size::Int  # Base population size multiplier
    memory_size::Int  # Memory size for adaptation
    
    # Population management
    n_inds_front::Int
    n_inds_front_max::Int
    
    # Adaptation parameters
    success_rate::Float64
    eb_hybrid_rate::Float64
    memory_iter::Int
    success_filled::Int
    
    # Memory arrays for adaptation
    memory_cr::Vector{Float64}
    memory_f::Vector{Float64}
    temp_success_cr::Vector{Float64}
    temp_success_f::Vector{Float64}
    fit_delta::Vector{Float64}
    weights::Vector{Float64}
    
    # Hybrid strategy arrays
    fit_mass::Vector{Float64}
    fit_temp::Vector{Float64}
    eb_hybrid_flag::Vector{Float64}
    
    # Temporary arrays
    trial::Vector{Float64}
    indices::Vector{Int}
    indices2::Vector{Int}
end

"""
    RDEx(;
        population_size::Int = 18,
        memory_size::Int = 5,
        information = Information(),
        options = Options()
    ) -> Algorithm

Create an RDEx algorithm instance compatible with Metaheuristics.jl.

# Arguments
- `population_size::Int = 18`: Base population size multiplier. Final population size will be `population_size * dimension`.
- `memory_size::Int = 5`: Memory size for parameter adaptation.

# Returns
An `Algorithm` object that can be used with `Metaheuristics.optimize`.

# Example
```julia
using RDEx
using Metaheuristics

f(x) = sum(x.^2)
bounds = [-5.0 -5.0; 5.0 5.0]  # 2D problem
algorithm = RDEx(population_size=18)
result = optimize(f, bounds, algorithm)
```
"""
function RDEx(;
    population_size::Int = 18,
    memory_size::Int = 5,
    information = Information(),
    options = Options()
)
    # Initialize with empty arrays - will be resized in initialize!
    parameters = RDEx(
        0,  # N
        population_size,  # population_size
        memory_size,  # memory_size
        0,  # n_inds_front
        0,  # n_inds_front_max
        0.5,  # success_rate
        Constants.EB_HYBRID_RATE_INIT,  # eb_hybrid_rate
        0,  # memory_iter
        0,  # success_filled
        Float64[],  # memory_cr
        Float64[],  # memory_f
        Float64[],  # temp_success_cr
        Float64[],  # temp_success_f
        Float64[],  # fit_delta
        Float64[],  # weights
        Float64[],  # fit_mass
        Float64[],  # fit_temp
        Float64[],  # eb_hybrid_flag
        Float64[],  # trial
        Int[],  # indices
        Int[]  # indices2
    )
    
    return Algorithm(
        parameters,
        information = information,
        options = options,
    )
end

"""
    initialize!(
        status,
        parameters::RDEx,
        problem,
        information,
        options,
        args...;
        kargs...
    ) -> State

Initialize the RDEx algorithm and create the initial Metaheuristics State.
"""
function initialize!(
    status,
    parameters::RDEx,
    problem,
    information,
    options,
    args...;
    kargs...
)
    D = getdim(problem)
    
    # Calculate population size
    initial_pop_size = parameters.population_size * D
    parameters.N = initial_pop_size
    parameters.n_inds_front = initial_pop_size
    parameters.n_inds_front_max = initial_pop_size
    
    # Set default options if needed
    if options.iterations == 0
        options.iterations = 500
    end
    
    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end
    
    # Initialize memory arrays
    parameters.memory_size = Constants.MEMORY_SIZE
    parameters.memory_iter = 0
    parameters.success_filled = 0
    parameters.success_rate = 0.5
    parameters.eb_hybrid_rate = Constants.EB_HYBRID_RATE_INIT
    
    # Allocate arrays
    popul_size = initial_pop_size * 2
    parameters.memory_cr = ones(parameters.memory_size)
    parameters.memory_f = ones(parameters.memory_size)
    parameters.temp_success_cr = zeros(popul_size)
    parameters.temp_success_f = zeros(popul_size)
    parameters.fit_delta = zeros(popul_size)
    parameters.weights = zeros(popul_size)
    parameters.fit_mass = zeros(parameters.n_inds_front_max)
    parameters.fit_temp = zeros(parameters.n_inds_front_max)
    parameters.eb_hybrid_flag = zeros(parameters.n_inds_front_max)
    parameters.trial = zeros(D)
    parameters.indices = zeros(Int, popul_size)
    parameters.indices2 = zeros(Int, popul_size)
    
    # Create initial State using gen_initial_state
    state = gen_initial_state(problem, parameters, information, options, status)
    
    # Sort population to ensure best is first
    if !isempty(state.population)
        sort!(state.population, lt=is_better)
        state.best_sol = state.population[1]
    end
    
    # Initialize fit_mass with current fitness values
    if !isempty(state.population)
        fits = fvals(state.population)
        n_front = min(length(fits), length(parameters.fit_mass))
        parameters.fit_mass[1:n_front] = fits[1:n_front]
    end
    
    state.iteration = 0
    state.start_time = time()
    
    return state
end

"""
    update_state!(
        status,
        parameters::RDEx,
        problem,
        information,
        options,
        args...;
        kargs...
    )

Update the state by running one generation of RDEx algorithm.
"""
function update_state!(
    status,
    parameters::RDEx,
    problem,
    information,
    options,
    args...;
    kargs...
)
    D = getdim(problem)
    rng = options.rng
    
    # Get bounds
    lower_bounds = problem.search_space.lb
    upper_bounds = problem.search_space.ub
    
    # Sort population by fitness
    sort!(status.population, lt=is_better)
    status.best_sol = status.population[1]
    
    # Extract positions and fitness values
    X = positions(status.population)
    fits = fvals(status.population)
    n_front = length(status.population)
    parameters.n_inds_front = n_front
    
    # Ensure arrays are properly sized
    if length(parameters.fit_temp) < n_front
        resize!(parameters.fit_temp, n_front)
    end
    if length(parameters.eb_hybrid_flag) < n_front
        resize!(parameters.eb_hybrid_flag, n_front)
    end
    if length(parameters.fit_mass) < n_front
        resize!(parameters.fit_mass, n_front)
    end
    
    # Prepare selection weights
    fit_temp2 = [exp(-i / n_front * 3) for i = 0:(n_front-1)]
    psizeval = max(2, Int(round(n_front * 0.7 * exp(-parameters.success_rate * 7))))
    fit_temp_prand = [3.0 * (n_front - i) for i = 0:(n_front-1)]
    psizeval2 = Int(round(n_front * 0.17 * (1 - 0.5 * status.f_calls / options.f_calls_limit)))
    if psizeval2 <= 1
        psizeval2 = 2
    end
    
    # Prepare indices for sorting
    if length(parameters.indices) < n_front
        resize!(parameters.indices, n_front)
        resize!(parameters.indices2, n_front)
    end
    perm = sortperm(fits[1:n_front])
    for i = 1:n_front
        parameters.indices[i] = perm[i] - 1  # 0-indexed
        parameters.indices2[i] = perm[i] - 1
    end
    
    # Calculate mean_f and sigma_f for standard mutation
    mean_f = 0.4 + tanh(parameters.success_rate * 5) * 0.25
    sigma_f = 0.02
    
    # Reset success tracking
    parameters.success_filled = 0
    max_success = length(parameters.temp_success_cr)
    
    # Process each individual in front
    new_solutions = Vector{typeof(status.best_sol)}()
    
    for ind_iter = 1:n_front
        # Check if we've exceeded evaluation limit
        if status.f_calls >= options.f_calls_limit
            break
        end
        
        # Select the chosen one (0-indexed)
        the_chosen_one = rand(rng, 0:(n_front-1))
        memory_current_index = rand(rng, 0:(parameters.memory_size-1))
        memory_current_index2 = rand(rng, 0:parameters.memory_size)
        
        # Decide on hybrid strategy
        rand_eb = rand(rng)
        if status.f_calls / options.f_calls_limit < 0.7
            rand_eb = 2.0
        end
        
        x_chosen = X[the_chosen_one+1, :]
        fit_chosen = fits[the_chosen_one+1]
        
        if (rand_eb * (1 - status.f_calls / options.f_calls_limit)) < parameters.eb_hybrid_rate
            # Elite-based mutation
            parameters.eb_hybrid_flag[ind_iter] = 1.0
            
            # Select indices for elite-based mutation
            prand_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand[1:psizeval2]) + 1
            prand = parameters.indices[prand_idx]
            while prand == the_chosen_one
                prand_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand[1:psizeval2]) + 1
                prand = parameters.indices[prand_idx]
            end
            
            rand1_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand) + 1
            rand1 = parameters.indices2[rand1_idx]
            while rand1 == prand
                rand1_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand) + 1
                rand1 = parameters.indices2[rand1_idx]
            end
            
            rand2_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand) + 1
            rand2 = parameters.indices[rand2_idx]
            while rand2 == prand || rand2 == rand1
                rand2_idx = SelectionUtils.weighted_sample(rng, fit_temp_prand) + 1
                rand2 = parameters.indices[rand2_idx]
            end
            
            # Get positions for ordering
            pos1 = X[prand+1, :]
            pos1_fit = fits[prand+1]
            pos3 = X[rand1+1, :]
            pos3_fit = fits[rand1+1]
            pos4 = X[rand2+1, :]
            pos4_fit = fits[rand2+1]
            
            # Order positions
            ordering = MutationUtils.elite_based_order(pos1, pos1_fit, pos3, pos3_fit, pos4, pos4_fit)
            
            # Generate F and Cr for elite-based
            if memory_current_index2 < parameters.memory_size
                F = rand(rng, Cauchy(parameters.memory_f[memory_current_index2+1], 0.1))
            else
                F = rand(rng, Cauchy(0.9, 0.1))
            end
            while F < 0.0
                if memory_current_index2 < parameters.memory_size
                    F = rand(rng, Cauchy(parameters.memory_f[memory_current_index2+1], 0.1))
                else
                    F = rand(rng, Cauchy(0.9, 0.1))
                end
            end
            if F > 1.0
                F = 1.0
            end
            
            if status.f_calls / options.f_calls_limit < 0.6 && F > 0.7
                F = 0.7
            end
            
            if memory_current_index2 < parameters.memory_size
                if parameters.memory_cr[memory_current_index2+1] < 0
                    Cr = 0.0
                else
                    Cr = rand(rng, Normal(parameters.memory_cr[memory_current_index2+1], 0.1))
                end
            else
                Cr = rand(rng, Normal(0.9, 0.1))
            end
            Cr = clamp(Cr, 0.0, 1.0)
            
            if status.f_calls / options.f_calls_limit < 0.25
                Cr = max(Cr, 0.7)
            end
            if status.f_calls / options.f_calls_limit < 0.5
                Cr = max(Cr, 0.6)
            end
            
            # Create trial
            trial, actual_cr = MutationUtils.create_elite_based_trial(
                x_chosen, ordering, F, Cr, lower_bounds, upper_bounds, rng
            )
        else
            # Standard mutation
            parameters.eb_hybrid_flag[ind_iter] = 0.0
            
            # Select indices for standard mutation
            prand_idx = rand(rng, 1:psizeval)
            prand = parameters.indices[prand_idx]
            while prand == the_chosen_one
                prand_idx = rand(rng, 1:psizeval)
                prand = parameters.indices[prand_idx]
            end
            
            rand1_idx = SelectionUtils.weighted_sample(rng, fit_temp2) + 1
            rand1 = parameters.indices2[rand1_idx]
            while rand1 == prand
                rand1_idx = SelectionUtils.weighted_sample(rng, fit_temp2) + 1
                rand1 = parameters.indices2[rand1_idx]
            end
            
            rand2_idx = rand(rng, 1:n_front)
            rand2 = parameters.indices[rand2_idx]
            while rand2 == prand || rand2 == rand1
                rand2_idx = rand(rng, 1:n_front)
                rand2 = parameters.indices[rand2_idx]
            end
            
            # Generate F and Cr for standard
            F = rand(rng, Normal(mean_f, sigma_f))
            while F < 0.0 || F > 1.0
                F = rand(rng, Normal(mean_f, sigma_f))
            end
            
            Cr = rand(rng, Normal(parameters.memory_cr[memory_current_index+1], 0.05))
            Cr = clamp(Cr, 0.0, 1.0)
            
            # Create trial
            trial, actual_cr = MutationUtils.create_standard_trial(
                x_chosen, X[prand+1, :], X[rand1+1, :], X[rand2+1, :],
                F, Cr, lower_bounds, upper_bounds, rng
            )
        end
        
        # Evaluate trial
        trial_sol = create_solution(trial, problem)
        trial_fit = trial_sol.f
        parameters.fit_temp[ind_iter] = trial_fit
        
        # Update if better
        if trial_fit <= fit_chosen
            status.population[the_chosen_one+1] = trial_sol
            if parameters.success_filled < max_success
                parameters.temp_success_cr[parameters.success_filled+1] = actual_cr
                parameters.temp_success_f[parameters.success_filled+1] = F
                parameters.fit_delta[parameters.success_filled+1] = abs(fit_chosen - trial_fit)
            end
            parameters.success_filled += 1
        end
    end
    
    # Update hybrid parameter
    fits_updated = fvals(status.population)
    AdaptationUtils.update_eb_hybrid_param!(
        parameters, parameters.eb_hybrid_flag, fits[1:n_front], 
        parameters.fit_temp, n_front
    )
    parameters.fit_mass[1:n_front] = fits_updated[1:n_front]
    
    # Update success rate and population size
    parameters.success_rate = parameters.success_filled / n_front
    new_n_inds_front = Int(round((4 - parameters.n_inds_front_max) / options.f_calls_limit * status.f_calls + parameters.n_inds_front_max))
    new_n_inds_front = clamp(new_n_inds_front, 4, parameters.n_inds_front_max)
    
    # Sort and truncate population if needed
    sort!(status.population, lt=is_better)
    if length(status.population) > new_n_inds_front
        status.population = status.population[1:new_n_inds_front]
    end
    parameters.n_inds_front = length(status.population)
    
    # Update memory
    AdaptationUtils.update_memory!(
        parameters, parameters.temp_success_cr, parameters.temp_success_f, parameters.fit_delta
    )
    
    # Update best solution
    status.best_sol = status.population[1]
    
    # Update success rate in status
    status.success_rate = parameters.success_rate
    
    return
end

"""
    final_stage!(
        status,
        parameters::RDEx,
        problem,
        information,
        options,
        args...;
        kargs...
    )

Finalize the optimization process.
"""
function final_stage!(
    status,
    parameters::RDEx,
    problem,
    information,
    options,
    args...;
    kargs...
)
    # Ensure best solution is correctly set
    if !isempty(status.population)
        sort!(status.population, lt=is_better)
        status.best_sol = status.population[1]
    end
    
    status.final_time = time()
    
    return
end

end # module MetaheuristicsAPI

