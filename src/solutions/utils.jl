_fix_type(x::AbstractVector{T}, bounds::BoxConstrainedSpace{T}) where T = x

function _fix_type(
        x::AbstractVector{<:Real},
        bounds::BoxConstrainedSpace{T},
    ) where T <:Integer
    round.(T, x)
end

function _fix_type(
        x::AbstractVector{<:Integer},
        bounds::BoxConstrainedSpace{T},
    ) where T <: AbstractFloat
    T.(x)
end

_fix_type(x, space) = x

