"""
Parameter adaptation mechanisms for RDEx optimizer.

Handles adaptation of F, Cr, and hybrid strategy parameters.
"""
module AdaptationUtils

import ..Constants

"""
    weighted_mean(values::Vector{Float64}, weights::Vector{Float64}, success_filled::Int, weights_storage::Vector{Float64}) -> Float64

Compute weighted mean of values using weights.
"""
function weighted_mean(
    values::Vector{Float64},
    weights::Vector{Float64},
    success_filled::Int,
    weights_storage::Vector{Float64}
)::Float64
    sum_weight = sum(weights[1:success_filled])
    if sum_weight == 0.0
        return 1.0
    end
    
    weights_storage[1:success_filled] = weights[1:success_filled] / sum_weight
    
    sum_square = sum(weights_storage[i] * values[i] * values[i] for i = 1:success_filled)
    sum_val = sum(weights_storage[i] * values[i] for i = 1:success_filled)
    
    if abs(sum_val) > 1e-8
        return sum_square / sum_val
    else
        return 1.0
    end
end

"""
    update_memory!(
        parameters,
        temp_success_cr::Vector{Float64},
        temp_success_f::Vector{Float64},
        fit_delta::Vector{Float64}
    )

Update memory for Cr and F parameters based on successful trials.
"""
function update_memory!(
    parameters,
    temp_success_cr::Vector{Float64},
    temp_success_f::Vector{Float64},
    fit_delta::Vector{Float64}
)::Nothing
    if parameters.success_filled != 0
        parameters.memory_cr[parameters.memory_iter+1] = 0.5 * (
            weighted_mean(temp_success_cr, fit_delta, parameters.success_filled, parameters.weights) + 
            parameters.memory_cr[parameters.memory_iter+1]
        )
        parameters.memory_f[parameters.memory_iter+1] = weighted_mean(
            temp_success_f, fit_delta, parameters.success_filled, parameters.weights
        )
        parameters.memory_iter = (parameters.memory_iter + 1) % parameters.memory_size
    end
    return nothing
end

"""
    update_eb_hybrid_param!(
        parameters,
        eb_hybrid_flag::Vector{Float64},
        fit_arr_front::Vector{Float64},
        fit_temp::Vector{Float64},
        n_inds_front::Int
    )

Update the elite-based hybrid rate parameter based on performance.
"""
function update_eb_hybrid_param!(
    parameters,
    eb_hybrid_flag::Vector{Float64},
    fit_arr_front::Vector{Float64},
    fit_temp::Vector{Float64},
    n_inds_front::Int
)::Nothing
    sum_eb_delta_fit = 0.0
    sum_origin_delta_fit = 0.0
    
    for chosen_one = 1:n_inds_front
        if eb_hybrid_flag[chosen_one] == 1.0
            if fit_temp[chosen_one] <= fit_arr_front[chosen_one]
                sum_eb_delta_fit += fit_arr_front[chosen_one] - fit_temp[chosen_one]
            end
        else
            if fit_temp[chosen_one] <= fit_arr_front[chosen_one]
                sum_origin_delta_fit += fit_arr_front[chosen_one] - fit_temp[chosen_one]
            end
        end
    end
    
    if sum_eb_delta_fit != 0.0 && sum_origin_delta_fit != 0.0
        parameters.eb_hybrid_rate = sum_eb_delta_fit / (sum_eb_delta_fit + sum_origin_delta_fit)
        parameters.eb_hybrid_rate = clamp(parameters.eb_hybrid_rate, 0.0, 1.0)
    else
        parameters.eb_hybrid_rate = Constants.EB_HYBRID_RATE_INIT
    end
    
    return nothing
end

end # module AdaptationUtils

