abstract type Neighborhood end
abstract type InternalNeighborhood <: Neighborhood end

struct NeighborhoodIterator{X, N} <: InternalNeighborhood
    x::X
    neighborhood::N
end

Base.@kwdef struct TwoOptNeighborhood <: Neighborhood
    k::Int = 2
end

function neighborhood_structure(x, s::Neighborhood, i)
    # The i-th neighbour in the k-th neighborhood around x
    n = nameof(typeof(s))

    println("""
            Define your neighborhood as follows:

            `
            function neighborhood_structure(x, s::$n, i)
                # ...
            end
            `
            """
           )
end


function neighborhood_structure(permutation, s::TwoOptNeighborhood, i)
    k = s.k-1
    if i isa Integer
        if !(1 <= i <= length(permutation)-k && 1 <= k <= length(permutation))
            return nothing
        end
    else
        i = rand(i, 1:length(permutation)-k)
    end
    reverse!(view(permutation, i:i+k))
    permutation
end


function Base.iterate(iter::NeighborhoodIterator, state=(copy(iter.x), 1))
    # reuse same x from state
    x, i = state
    x .= iter.x
    v = neighborhood_structure(x, iter.neighborhood, i)
    if isnothing(v)
        return v
    end
    # return to current solution
    v, (x, i+1)
end

