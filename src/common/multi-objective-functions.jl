"""
    gen_ref_dirs(dimension, n_paritions)
Generates Das and Dennis's structured reference points. `dimension` could be
the number of objective functions in multi-objective functions.
"""
function gen_ref_dirs(dimension, n_paritions)
    return  gen_weights(dimension, n_paritions)
end

function gen_weights(a, b)
    nobj = a;
    H    = b;
    a    = zeros(nobj);
    d    = H;
    w    = [];
    produce_weight!(a, 1, d, H, nobj, w)
    return Array.(w)
end

function  produce_weight!(a, i, d, H, nobj, w)
    for k=0:d
        if i<nobj
            a[i] = k;
            d2   = d - k;
            produce_weight!(a, i+1, d2, H, nobj, w);
        else
            a[i] = d;
            push!(w, a/H)
            break;
        end
    end
end


"""
    ideal(points)

Computes the ideal point from a provided array of `Vector`s or a population or row vectors
in a `Matrix`.
"""
function ideal(points::Array{Vector{T}}) where T <: Real 

    (isempty(points) || isempty(points[1])) && return zeros(0)
    
    ideal = points[1]

    for point in points
        ideal = min.(ideal, point)
    end

    return ideal

end


"""
    nadir(points)

Computes the nadir point from a provided array of `Vector`s or a population or row vectors
in a `Matrix`.
"""
function nadir(points::Array{Vector{T}})  where T <: Real 
    (isempty(points) || isempty(points[1])) && return zeros(0)
    nadir = points[1]

    for point in points
        nadir = max.(nadir, point)
    end

    return nadir
end


function ideal(population::Array{T}) where T <: xFgh_solution
    isempty(population) && (return zeros(0))

    mask = sum_violations.(population) .== 0
    if count(mask) == 0
        @warn "Ideal point was computed using infeasible solutions. Use `ideal(fvals(population))` to ignore feasibility."
        return ideal(fval.(population))
    end
    

    ideal(fval.(population[mask]))
end
ideal(A::Matrix) = ideal([A[i,:]  for i in 1:size(A,1)])

function nadir(population::Array{T}) where T <: xFgh_solution
    isempty(population) && (return zeros(0))
    mask = sum_violations.(population) .== 0

    if count(mask) == 0
        @warn "Nadir point was computed using infeasible solutions. Use `nadir(fvals(population))` to ignore feasibility."
        return nadir(fval.(population))
    end

    nadir(fval.(population[mask]))
end

nadir(A::Matrix) = nadir([A[i,:]  for i in 1:size(A,1)])
