abstract type Neighbourhood end
abstract type InternalNeighbourhood <: Neighbourhood end

struct NeighborhoodIterator{X, N} <: InternalNeighbourhood
    x::X
    neighborhood::N
end

Base.@kwdef struct TwoOpt
    x::Vector{Int}
    k::Int = 2
end

function Base.iterate(two_opt::TwoOpt, state=1)
    i, permutation, k = state, two_opt.x, two_opt.k-1
    if !(1 <= i <= length(permutation)-k && 1 <= k <= length(permutation))
        return nothing
    end
    reverse!(view(permutation, i:i+k))
    permutation, state+1
end


function neighborhood_structure(x, s::Neighbourhood, i)
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


#=
function local_search(x, neighborhood::TwoOpt, problem)
    # FIXME: overwriting best.x by every xnew
    best = MH.create_solution(x, problem)
    for xnew in neighborhood
        sol = MH.create_solution(xnew, problem)
        if isnothing(best) || MH.is_better(sol, best)
            best = sol
        end
    end
    best
end
=#
