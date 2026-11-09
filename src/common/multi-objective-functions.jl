"""
    gen_ref_dirs(dimension, n_paritions)
Generates Das and Dennis's structured reference points. `dimension` could be
the number of objective functions in multi-objective functions.
"""
function gen_ref_dirs(dimension, n_paritions)
    return gen_weights(dimension, n_paritions)
end

function gen_weights(a, b)
    nobj = a
    H = b
    a = zeros(nobj)
    d = H
    w = Vector{Float64}[]
    produce_weight!(a, 1, d, H, nobj, w)
    return w
end

function produce_weight!(a, i, d, H, nobj, w)
    for k in 0:d
        if i < nobj
            a[i] = k
            d2 = d - k
            produce_weight!(a, i + 1, d2, H, nobj, w)
        else
            a[i] = d
            push!(w, a / H)
            break
        end
    end
end

"""
    ideal(points)

Computes the ideal point from a provided array of `Vector`s or a population or row vectors
in a `Matrix`.
"""
function ideal(points::Array{Vector{T}}) where {T <: Real}

    (isempty(points) || isempty(points[1])) && return zeros(0)

    ideal = copy(points[1])

    for point in points
        ideal .= min.(ideal, point)
    end

    return ideal

end

"""
    nadir(points)

Computes the nadir point from a provided array of `Vector`s or a population or row vectors
in a `Matrix`.
"""
function nadir(points::Array{Vector{T}}) where {T <: Real}
    (isempty(points) || isempty(points[1])) && return zeros(0)
    nadir = copy(points[1])

    for point in points
        nadir .= max.(nadir, point)
    end

    return nadir
end

function ideal(population::Array{T}) where {T <: xFgh_solution}
    isempty(population) && (return zeros(0))

    return _extreme_point(min, population, "Ideal", "ideal")
end
ideal(A::Matrix) = _extreme_point(min, A)

function nadir(population::Array{T}) where {T <: xFgh_solution}
    isempty(population) && (return zeros(0))
    return _extreme_point(max, population, "Nadir", "nadir")
end

nadir(A::Matrix) = _extreme_point(max, A)

# component-wise min/max over feasible solutions (all if none is feasible)
function _extreme_point(op::F, population::Array{<:xFgh_solution}, name, fname) where {F}
    p = nothing
    for s in population
        sum_violations(s) == 0 || continue
        p = isnothing(p) ? copy(fval(s)) : (p .= op.(p, fval(s)))
    end
    isnothing(p) || return p

    @warn "$name point was computed using infeasible solutions. Use `$fname(fvals(population))` to ignore feasibility."
    p = copy(fval(first(population)))
    for s in population
        p .= op.(p, fval(s))
    end
    return p
end

function _extreme_point(op::F, A::Matrix) where {F}
    (size(A, 1) == 0 || size(A, 2) == 0) && return zeros(0)
    p = A[1, :]
    for j in axes(A, 2), i in axes(A, 1)
        @inbounds p[j] = op(p[j], A[i, j])
    end
    return p
end
