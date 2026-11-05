"""
Type definitions for RDEx optimizer.
"""
module Types

import ..RandomUtils: RNGState

"""
    Optimizer

Main optimizer state structure for RDEx algorithm.

All state is encapsulated here, avoiding global variables.
"""
mutable struct Optimizer
    # Objective function and problem definition
    objective_func::Function
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    dimension::Int
    
    # Evaluation tracking
    n_evaluations::Int
    max_evaluations::Int
    
    # Random number generation
    rng::RNGState
    
    # Memory parameters for adaptation
    memory_size::Int
    memory_iter::Int
    success_filled::Int
    memory_current_index::Int
    memory_current_index2::Int
    
    # Population management
    n_inds_current::Int
    n_inds_front::Int
    n_inds_front_max::Int
    new_n_inds_front::Int
    popul_size::Int
    
    # Selection indices
    the_chosen_one::Int
    pf_index::Int
    
    # Best solution tracking
    best_fitness::Float64
    best_solution::Vector{Float64}
    
    # Adaptation parameters
    success_rate::Float64
    eb_hybrid_rate::Float64
    
    # Population arrays
    popul::Matrix{Float64}
    popul_front::Matrix{Float64}
    popul_temp::Matrix{Float64}
    fit_arr::Vector{Float64}
    fit_arr_copy::Vector{Float64}
    fit_arr_front::Vector{Float64}
    trial::Vector{Float64}
    
    # Memory and adaptation
    temp_success_cr::Vector{Float64}
    temp_success_f::Vector{Float64}
    memory_cr::Vector{Float64}
    memory_f::Vector{Float64}
    fit_delta::Vector{Float64}
    weights::Vector{Float64}
    
    # Indices for sorting
    indices::Vector{Int}
    indices2::Vector{Int}
    
    # Hybrid strategy
    fit_mass::Vector{Float64}
    fit_temp::Vector{Float64}
    eb_hybrid_flag::Vector{Float64}
end

"""
    OptimizationResult

Result structure returned by the optimize function.
"""
struct OptimizationResult
    best_solution::Vector{Float64}
    best_fitness::Float64
    n_evaluations::Int
    convergence_history::Vector{Float64}
end

end # module Types

