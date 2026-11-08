function reflected_back_to_bound!(x, bounds::BoxConstrainedSpace)
    !bounds.rigid && (return x)

    for i in 1:getdim(bounds)
        l = bounds.lb[i]
        u = bounds.ub[i]
        xi = x[i]
        l <= xi <= u && continue
        # one reflection (most common case)
        y = xi < l ? 2l - xi : 2u - xi
        x[i] = l <= y <= u ? y : _reflect_far(xi, l, u)
    end
    x
end

# closed form of repeated reflections (triangle wave of period 2(u - l))
function _reflect_far(xi, l, u)
    ρ = u - l
    (isnan(xi) || iszero(ρ)) && return oftype(xi, l)
    isinf(xi) && return oftype(xi, xi < l ? l : u)
    t = mod(xi - l, 2ρ)
    return clamp(t <= ρ ? l + t : l + (2ρ - t), l, u)
end

reflected_back_to_bound!(x, bounds::AbstractMatrix) = reflected_back_to_bound!(
    x,
    _mat_to_bounds(bounds),
)

function replace_with_random_in_bounds!(
    x,
    bounds::BoxConstrainedSpace,
    rng = default_rng_mh(),
)
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    Δ = bounds.Δ
    for i in 1:getdim(bounds)
        if !(lb[i] <= x[i] <= ub[i])
            x[i] = lb[i] + rand(rng) * Δ[i]
        end
    end
    x
end

replace_with_random_in_bounds!(x, bounds::AbstractMatrix) = replace_with_random_in_bounds!(
    x,
    _mat_to_bounds(bounds),
)

function wrap_to_bounds!(x, bounds::BoxConstrainedSpace)
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    for i in 1:getdim(bounds)
        l = lb[i]
        u = ub[i]
        if !(l <= x[i] <= u)
            ρ = u - l
            x[i] = x[i] < l ? u - (l - x[i]) % ρ : l + (x[i] - u) % ρ
        end

    end
    x
end
wrap_to_bounds!(x, bounds::AbstractMatrix) = wrap_to_bounds!(x, _mat_to_bounds(bounds))

function reset_to_violated_bounds!(x, bounds::BoxConstrainedSpace)
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
reset_to_violated_bounds!(x, bounds::AbstractMatrix) = reset_to_violated_bounds!(
    x,
    _mat_to_bounds(bounds),
)

function evo_boundary_repairer!(
    x,
    x_best,
    bounds::BoxConstrainedSpace,
    rng = default_rng_mh(),
)
    !bounds.rigid && (return x)
    lb = bounds.lb
    ub = bounds.ub
    for i in 1:getdim(bounds)
        l = lb[i]
        u = ub[i]
        if l > x[i]
            α = rand(rng)
            x[i] = α * l + (1.0 - α) * x_best[i]
        elseif x[i] > u
            β = rand(rng)
            x[i] = β * u + (1.0 - β) * x_best[i]
        end

    end
    x
end
evo_boundary_repairer!(
    x,
    x_best,
    bounds::AbstractMatrix,
    rng = default_rng_mh(),
) = evo_boundary_repairer!(x, x_best, _mat_to_bounds(bounds), rng)

function is_in_bounds(x, bounds::BoxConstrainedSpace)
    for i in 1:getdim(bounds)
        l = bounds.lb[i]
        u = bounds.ub[i]
        if !(l <= x[i] <= u)
            return false
        end

    end
    return true
end

is_in_bounds(x, bounds::AbstractMatrix) = is_in_bounds(x, _mat_to_bounds(bounds))
