function reflected_back_to_bound!(x, bounds)
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        while !( l <= x[i] <= u )
            x[i] = x[i] < l ? 2l - x[i] : 2u - x[i]
        end

    end
    x
end

function replace_with_random_in_bounds!(x, bounds) 
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        if !( l <= x[i] <= u )
            x[i] = l + rand() * (u - l)
        end

    end
    x
end

function wrap_to_bounds!(x, bounds)
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        if !( l <= x[i] <= u )
            ρ = u - l
            x[i] = x[i] < l ? u - (l - x[i]) % ρ : l + (x[i] - u) % ρ
        end

    end
    x
end

function reset_to_violated_bounds!(x, bounds)
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        if l > x[i]
            x[i] = l
        elseif x[i] > u 
            x[i] = u
        end

    end
    x 
end

function evo_boundary_repairer!(x, x_best, bounds)            
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        if l > x[i]
            α = rand()
            x[i] =  α*l + (1.0 - α) * x_best[i]
        elseif x[i] > u  
            β = rand()
            x[i] =  β*u + (1.0 - β) * x_best[i]
        end

    end
    x  
end


function is_in_bounds(x, bounds) 
    for i in 1:size(bounds,2)
        l = bounds[1, i]
        u = bounds[2, i]
        if !( l <= x[i] <= u )
            return false
        end

    end
    return true
end

