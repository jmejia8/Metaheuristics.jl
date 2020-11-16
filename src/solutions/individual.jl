mutable struct xf_indiv <: AbstractSolution # Single Objective
    x::Vector{Float64}
    f::Float64
end

mutable struct xfg_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
end

mutable struct xfgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Float64
    g::Vector{Float64}
    h::Vector{Float64}
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|

end

function xfgh_indiv(
    x::Vector{Float64},
    f::Float64,
    g::Vector{Float64},
    h::Vector{Float64};
    sum_violations = 0
)
    if sum_violations <= 0
        sum_violations = violationsSum(g, h)
    end

    xfgh_indiv(x, f, g, h, sum_violations)
end

mutable struct xFgh_indiv # Single Objective Constraied
    x::Vector{Float64}
    f::Vector{Float64}
    g::Vector{Float64}
    h::Vector{Float64}
    rank::Int
    crowding::Float64
    sum_violations::Float64 # ∑ max(0,g) + ∑|h|
end

function xFgh_indiv(
    x::Vector{Float64},
    f::Vector{Float64},
    g::Vector{Float64},
    h::Vector{Float64};
    rank = 0,
    crowding = 0.0,
    sum_violations = 0.0
)

    if sum_violations <= 0
        sum_violations = violationsSum(g, h)
    end
    xFgh_indiv(x, f, g, h, Int(rank), crowding, sum_violations)
end



############################################################
# Generate solutions depending on the objective function
# output
#############################################################

function generateChild(x::Vector{Float64}, fResult::Float64)
    return xf_indiv(x, fResult)
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Float64,Array{Float64,1}})
    f, g = fResult
    return xfgh_indiv(x, f, g, [0.0])
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    f, g, h = fResult
    return xfgh_indiv(x, f, g, h)
end

function generateChild(x::Vector{Float64}, fResult::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}})
    f, g, h = fResult
    return xFgh_indiv(x, f, g, h)
end

function inferType(fVal::Tuple{Float64})
    return xf_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1}})
    return xfg_indiv
end

function inferType(fVal::Tuple{Float64,Array{Float64,1},Array{Float64,1}})
    return xfgh_indiv
end

function inferType(fVal::Tuple{Array{Float64,1},Array{Float64,1},Array{Float64,1}})
    return xFgh_indiv
end


# getters for the above structures

"""
    get_position(solution)

Get the position vector.
"""
get_position(solution::Union{xf_indiv, xfg_indiv, xFgh_indiv}) = solution.x

"""
    fval(solution)

Get the objective function value (fitness) of a solution.
"""
fval(solution::Union{xf_indiv, xfg_indiv, xfgh_indiv, xFgh_indiv}) = solution.f
