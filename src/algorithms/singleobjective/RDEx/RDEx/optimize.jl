"""
Main optimization cycle for RDEx algorithm.
"""
module OptimizeUtils

import ..Types
import ..RandomUtils
import ..SelectionUtils
import ..MutationUtils
import ..AdaptationUtils
import ..OptimizerUtils
import ..Constants

"""
    evaluate_initial_population!(opt::Types.Optimizer)

Evaluate the initial population and find the best solution.
"""
function evaluate_initial_population!(opt::Types.Optimizer)::Nothing
    opt.best_fitness = Inf
    for ind_iter = 1:opt.n_inds_front
        opt.fit_arr[ind_iter] = opt.objective_func(opt.popul[ind_iter, :])
        opt.n_evaluations += 1
        AdaptationUtils.find_and_save_best!(opt, ind_iter == 1, ind_iter - 1)
    end
    return nothing
end

"""
    prepare_selection_weights(opt::Types.Optimizer) -> Tuple{Vector{Float64}, Int, Vector{Float64}, Int}

Prepare selection weights and parameters for the current generation.
Returns (fit_temp2, psizeval, fit_temp_prand, psizeval2).
"""
function prepare_selection_weights(opt::Types.Optimizer)::Tuple{Vector{Float64}, Int, Vector{Float64}, Int}
    fit_temp2 = [exp(-i / opt.n_inds_front * 3) for i = 0:(opt.n_inds_front-1)]
    psizeval = max(2, Int(round(opt.n_inds_front * 0.7 * exp(-opt.success_rate * 7))))
    fit_temp_prand = [3.0 * (opt.n_inds_front - i) for i = 0:(opt.n_inds_front-1)]
    psizeval2 = Int(round(opt.n_inds_front * 0.17 * (1 - 0.5 * opt.n_evaluations / opt.max_evaluations)))
    if psizeval2 <= 1
        psizeval2 = 2
    end
    return (fit_temp2, psizeval, fit_temp_prand, psizeval2)
end

"""
    select_standard_indices(opt::Types.Optimizer, psizeval::Int, fit_temp2::Vector{Float64}) -> Tuple{Int, Int, Int}

Select indices for standard mutation strategy.
Returns (prand, rand1, rand2).
"""
function select_standard_indices(
    opt::Types.Optimizer,
    psizeval::Int,
    fit_temp2::Vector{Float64}
)::Tuple{Int, Int, Int}
    # Select prand
    prand = opt.indices[RandomUtils.int_random(opt.rng, psizeval)+1]
    while prand == opt.the_chosen_one
        prand = opt.indices[RandomUtils.int_random(opt.rng, psizeval)+1]
    end
    
    # Select Rand1
    rand1 = opt.indices2[SelectionUtils.weighted_sample(opt.rng, fit_temp2) + 1]
    while rand1 == prand
        rand1 = opt.indices2[SelectionUtils.weighted_sample(opt.rng, fit_temp2) + 1]
    end
    
    # Select Rand2
    rand2 = opt.indices[RandomUtils.int_random(opt.rng, opt.n_inds_front)+1]
    while rand2 == prand || rand2 == rand1
        rand2 = opt.indices[RandomUtils.int_random(opt.rng, opt.n_inds_front)+1]
    end
    
    return (prand, rand1, rand2)
end

"""
    select_elite_indices(opt::Types.Optimizer, fit_temp_prand::Vector{Float64}, psizeval2::Int) -> Tuple{Int, Int, Int}

Select indices for elite-based mutation strategy.
Returns (prand, rand1, rand2).
"""
function select_elite_indices(
    opt::Types.Optimizer,
    fit_temp_prand::Vector{Float64},
    psizeval2::Int
)::Tuple{Int, Int, Int}
    prand = opt.indices[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand[1:psizeval2]) + 1]
    while prand == opt.the_chosen_one
        prand = opt.indices[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand[1:psizeval2]) + 1]
    end
    
    rand1 = opt.indices2[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand) + 1]
    while rand1 == prand
        rand1 = opt.indices2[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand) + 1]
    end
    
    rand2 = opt.indices[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand) + 1]
    while rand2 == prand || rand2 == rand1
        rand2 = opt.indices[SelectionUtils.weighted_sample(opt.rng, fit_temp_prand) + 1]
    end
    
    return (prand, rand1, rand2)
end

"""
    generate_standard_F_Cr(opt::Types.Optimizer, mean_f::Float64, sigma_f::Float64) -> Tuple{Float64, Float64}

Generate F and Cr parameters for standard mutation.
Returns (F, Cr).
"""
function generate_standard_F_Cr(
    opt::Types.Optimizer,
    mean_f::Float64,
    sigma_f::Float64
)::Tuple{Float64, Float64}
    F = RandomUtils.norm_rand(opt.rng, mean_f, sigma_f)
    while F < 0.0 || F > 1.0
        F = RandomUtils.norm_rand(opt.rng, mean_f, sigma_f)
    end
    
    Cr = RandomUtils.norm_rand(opt.rng, opt.memory_cr[opt.memory_current_index+1], 0.05)
    Cr = clamp(Cr, 0.0, 1.0)
    
    return (F, Cr)
end

"""
    generate_elite_F_Cr(opt::Types.Optimizer) -> Tuple{Float64, Float64}

Generate F and Cr parameters for elite-based mutation.
Returns (F, Cr).
"""
function generate_elite_F_Cr(opt::Types.Optimizer)::Tuple{Float64, Float64}
    # Generate F for elite-based
    if opt.memory_current_index2 < opt.memory_size
        F = RandomUtils.cachy_rand(opt.rng, opt.memory_f[opt.memory_current_index2+1], 0.1)
    else
        F = RandomUtils.cachy_rand(opt.rng, 0.9, 0.1)
    end
    while F < 0.0
        if opt.memory_current_index2 < opt.memory_size
            F = RandomUtils.cachy_rand(opt.rng, opt.memory_f[opt.memory_current_index2+1], 0.1)
        else
            F = RandomUtils.cachy_rand(opt.rng, 0.9, 0.1)
        end
    end
    if F > 1.0
        F = 1.0
    end
    
    if opt.n_evaluations / opt.max_evaluations < 0.6 && F > 0.7
        F = 0.7
    end
    
    # Generate Cr for elite-based
    if opt.memory_current_index2 < opt.memory_size
        if opt.memory_cr[opt.memory_current_index2+1] < 0
            Cr = 0.0
        else
            Cr = RandomUtils.norm_rand(opt.rng, opt.memory_cr[opt.memory_current_index2+1], 0.1)
        end
    else
        Cr = RandomUtils.norm_rand(opt.rng, 0.9, 0.1)
    end
    Cr = clamp(Cr, 0.0, 1.0)
    
    if opt.n_evaluations / opt.max_evaluations < 0.25
        Cr = max(Cr, 0.7)
    end
    if opt.n_evaluations / opt.max_evaluations < 0.5
        Cr = max(Cr, 0.6)
    end
    
    return (F, Cr)
end

"""
    process_individual!(opt::Types.Optimizer, mean_f::Float64, sigma_f::Float64, 
                       fit_temp2::Vector{Float64}, psizeval::Int, 
                       fit_temp_prand::Vector{Float64}, psizeval2::Int)

Process a single individual in the current generation.
"""
function process_individual!(
    opt::Types.Optimizer,
    mean_f::Float64,
    sigma_f::Float64,
    fit_temp2::Vector{Float64},
    psizeval::Int,
    fit_temp_prand::Vector{Float64},
    psizeval2::Int
)::Nothing
    opt.the_chosen_one = RandomUtils.int_random(opt.rng, opt.n_inds_front)
    opt.memory_current_index = RandomUtils.int_random(opt.rng, opt.memory_size)
    opt.memory_current_index2 = RandomUtils.int_random(opt.rng, opt.memory_size + 1)
    
    # Decide on hybrid strategy
    rand_eb = RandomUtils.random_val(opt.rng, 0.0, 1.0)
    if opt.n_evaluations / opt.max_evaluations < 0.7
        rand_eb = 2.0
    end
    
    if (rand_eb * (1 - opt.n_evaluations / opt.max_evaluations)) < opt.eb_hybrid_rate
        # Elite-based mutation
        opt.eb_hybrid_flag[opt.the_chosen_one+1] = 1.0
        
        prand, rand1, rand2 = select_elite_indices(opt, fit_temp_prand, psizeval2)
        arch_ordering = MutationUtils.elite_based_order!(opt, prand, rand1, rand2)
        F, Cr = generate_elite_F_Cr(opt)
        actual_cr = MutationUtils.create_elite_based_trial!(opt, arch_ordering, F, Cr)
    else
        # Standard mutation
        opt.eb_hybrid_flag[opt.the_chosen_one+1] = 0.0
        
        prand, rand1, rand2 = select_standard_indices(opt, psizeval, fit_temp2)
        F, Cr = generate_standard_F_Cr(opt, mean_f, sigma_f)
        actual_cr = MutationUtils.create_standard_trial!(opt, prand, rand1, rand2, F, Cr)
    end
    
    # Evaluate trial
    temp_fit = opt.objective_func(opt.trial)
    opt.n_evaluations += 1
    opt.fit_temp[opt.the_chosen_one+1] = temp_fit
    
    # Update if better
    if temp_fit <= opt.fit_arr_front[opt.the_chosen_one+1]
        opt.popul[opt.n_inds_current + opt.success_filled + 1, :] = opt.trial
        opt.popul_front[opt.pf_index+1, :] = opt.trial
        opt.fit_arr[opt.n_inds_current + opt.success_filled + 1] = temp_fit
        opt.fit_arr_front[opt.pf_index+1] = temp_fit
        AdaptationUtils.find_and_save_best!(opt, false, opt.n_inds_current + opt.success_filled)
        opt.temp_success_cr[opt.success_filled+1] = actual_cr
        opt.temp_success_f[opt.success_filled+1] = F
        opt.fit_delta[opt.success_filled+1] = abs(opt.fit_arr_front[opt.the_chosen_one+1] - temp_fit)
        opt.success_filled += 1
        opt.pf_index = (opt.pf_index + 1) % opt.n_inds_front
    end
    
    return nothing
end

"""
    update_generation!(opt::Types.Optimizer)

Update parameters and population at the end of a generation.
"""
function update_generation!(opt::Types.Optimizer)::Nothing
    # Update hybrid parameter
    fit_temp_prand = opt.fit_temp[1:opt.n_inds_front]
    AdaptationUtils.update_eb_hybrid_param!(opt, opt.eb_hybrid_flag, opt.fit_mass, fit_temp_prand)
    opt.fit_mass[1:opt.n_inds_front] = opt.fit_arr_front[1:opt.n_inds_front]
    
    # Update success rate and population size
    opt.success_rate = opt.success_filled / opt.n_inds_front
    opt.new_n_inds_front = Int(round((4 - opt.n_inds_front_max) / opt.max_evaluations * opt.n_evaluations + opt.n_inds_front_max))
    AdaptationUtils.remove_worst!(opt, opt.n_inds_front, opt.new_n_inds_front)
    opt.n_inds_front = opt.new_n_inds_front
    AdaptationUtils.update_memory!(opt)
    opt.n_inds_current = opt.n_inds_front + opt.success_filled
    opt.success_filled = 0
    
    # Sort and truncate if needed
    OptimizerUtils.truncate_population!(opt)
    return nothing
end

"""
    optimize_cycle!(opt::Types.Optimizer)

Run the main optimization cycle until max_evaluations is reached.
"""
function optimize_cycle!(opt::Types.Optimizer)::Nothing
    # Initial evaluation
    evaluate_initial_population!(opt)
    
    # Sort initial population
    OptimizerUtils.sort_population!(opt, false)
    OptimizerUtils.initialize_front_population!(opt)
    
    # Main optimization loop
    while opt.n_evaluations < opt.max_evaluations
        mean_f = 0.4 + tanh(opt.success_rate * 5) * 0.25
        sigma_f = 0.02
        
        # Sort current and front populations
        OptimizerUtils.sort_population!(opt, false)
        OptimizerUtils.sort_population!(opt, true)
        
        # Prepare selection weights
        fit_temp2, psizeval, fit_temp_prand, psizeval2 = prepare_selection_weights(opt)
        
        # Process each individual in front
        for ind_iter = 1:opt.n_inds_front
            process_individual!(opt, mean_f, sigma_f, fit_temp2, psizeval, fit_temp_prand, psizeval2)
        end
        
        # Update generation
        update_generation!(opt)
    end
    
    return nothing
end

end # module OptimizeUtils

