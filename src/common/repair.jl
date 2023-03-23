_mat_to_bounds(bounds::AbstractMatrix) = Bounds(lb=bounds[1,:], ub=bounds[2,:])


function reflected_back_to_bound!(x, bounds::Bounds)
    !bounds.rigid && (return x)

    for i in 1:getdim(bounds)
        l = bounds.lb[i]
        u = bounds.ub[i]
        while !( l <= x[i] <= u )
            x[i] = x[i] < l ? 2l - x[i] : 2u - x[i]
        end

    end
    x
end

reflected_back_to_bound!(x, bounds::AbstractMatrix) = reflected_back_to_bound!(x, _mat_to_bounds(bounds))


function replace_with_random_in_bounds!(x, bounds::Bounds, rng = default_rng_mh())
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    Δ = bounds.Δ
    for i in 1:getdim(bounds)
        if !( lb[i] <= x[i] <= ub[i] )
            x[i] = lb[i] + rand(rng) * Δ[i]
        end
    end
    x
end

replace_with_random_in_bounds!(x, bounds::AbstractMatrix) = replace_with_random_in_bounds!(x, _mat_to_bounds(bounds))

function wrap_to_bounds!(x, bounds::Bounds)
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    for i in 1:getdim(bounds)
        l = lb[i]
        u = ub[i]
        if !( l <= x[i] <= u )
            ρ = u - l
            x[i] = x[i] < l ? u - (l - x[i]) % ρ : l + (x[i] - u) % ρ
        end

    end
    x
end
wrap_to_bounds!(x, bounds::AbstractMatrix) = wrap_to_bounds!(x, _mat_to_bounds(bounds))

function reset_to_violated_bounds!(x, bounds::Bounds)
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    for i in 1:getdim(bounds)
        l = lb[i]
        u = ub[i]
        if l > x[i]
            x[i] = l
        elseif x[i] > u 
            x[i] = u
        end

    end
    x 
end
reset_to_violated_bounds!(x, bounds::AbstractMatrix) = reset_to_violated_bounds!(x, _mat_to_bounds(bounds))

function evo_boundary_repairer!(x, x_best, bounds::Bounds, rng=default_rng_mh())
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    for i in 1:getdim(bounds)
        l = lb[i]
        u = ub[i]
        if l > x[i]
            α = rand(rng)
            x[i] =  α*l + (1.0 - α) * x_best[i]
        elseif x[i] > u  
            β = rand(rng)
            x[i] =  β*u + (1.0 - β) * x_best[i]
        end

    end
    x  
end
evo_boundary_repairer!(x, x_best, bounds::AbstractMatrix, rng=default_rng_mh()) = evo_boundary_repairer!(x, x_best, _mat_to_bounds(bounds), rng)


function is_in_bounds(x, bounds::Bounds) 
    for i in 1:getdim(bounds)
        l = bounds.lb[i]
        u = bounds.ub[i]
        if !( l <= x[i] <= u )
            return false
        end

    end
    return true
end

is_in_bounds(x, bounds::AbstractMatrix) = is_in_bounds(x, _mat_to_bounds(bounds))
