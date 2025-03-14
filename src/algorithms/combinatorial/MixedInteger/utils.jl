function _mixed_to_continous_space(mixed_space::MixedSpace)
    dim = getdim(mixed_space)
    lower_bounds = []
    upper_bounds = []

    rigid = true 
    for k in mixed_space.key_order
        space = mixed_space.domain[k]
        if space isa BoxConstrainedSpace
            append!(lower_bounds, space.lb)
            append!(upper_bounds, space.ub)
            if !space.rigid
                rigid = false
            end
        elseif space isa BitArraySpace
            append!(lower_bounds, zeros(Int, getdim(space)))
            append!(upper_bounds, ones( Int, getdim(space)))
        else
            error("Only BoxConstrainedSpace and BitArraySpace are currently supported.")
        end
    end
    lb = [v for v in promote(lower_bounds...)]
    ub = [v for v in promote(upper_bounds...)]
    BoxConstrainedSpace(lb, ub;rigid)
end

function _bound_type!(
        x::AbstractVector{<:Real},
        bounds::BoxConstrainedSpace{T},
    ) where T 
    # TODO: update following bounds due to issue "InexactError: trunc(Int64, <big value>)"
    clamp!(x, typemin(T)/2, typemax(T)/2)
    _fix_type(x, bounds)
end

function _promote_vals_and_types!(v, a, b, domain)
    x = v[a:b]
    if domain isa BoxConstrainedSpace
        # convert x into appropriate values
        x = _bound_type!(x, domain)
        reset_to_violated_bounds!(x, domain)
    elseif domain isa BitArraySpace
        # convert x into boolean
        x = clamp.(x, 0, 1) .> 0.5
    end
    v[a:b] = x
    x
end


function vec_to_dict(v::AbstractVector, mixed_space::MixedSpace)
    @assert length(v) == getdim(mixed_space)
    keys = mixed_space.key_order
    domain = mixed_space.domain

    last = cumsum(getdim(domain[k]) for k in keys)

    # create dictionaries
    decisions = Dict(k => begin
                            a, b = (last[i] .- getdim(domain[k])+1), last[i]
                            _promote_vals_and_types!(v, a, b, domain[k])
                     end
                     for (i, k) in enumerate(keys))

    decisions
end
