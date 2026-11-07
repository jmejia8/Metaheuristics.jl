"""
    OrderCrossover()
Order crossover for representation where order is important. Suitable for permutation
representation.

Each row must be a permutation of the same set of values; otherwise an `ArgumentError`
is thrown. If the number of rows is odd, the last row is copied unchanged.
"""
struct OrderCrossover{R}
    rng::R
end

OrderCrossover(; rng = default_rng_mh()) = OrderCrossover(rng)

function crossover(population, parameters::OrderCrossover)
    O = positions(population)
    N, D = size(O)
    s = rand(parameters.rng, 1:D, N) # slash points
    PA = Vector{eltype(O)}(undef, D) # parent A
    PB = Vector{eltype(O)}(undef, D) # parent B
    seen = falses(D)
    for i in 1:2:N - 1
        @inbounds for j in 1:D
            PA[j] = O[i, j]
            PB[j] = O[i + 1, j]
        end
        if _is_permutation_1toD!(seen, PA) && _is_permutation_1toD!(seen, PB)
            _order_fill!(O, i, PA, PB, s[i], seen)
            _order_fill!(O, i + 1, PB, PA, s[i], seen)
        else
            # values are not 1:D, use the generic (slower) version
            _order_fill_generic!(O, i, PA, PB, s[i])
            _order_fill_generic!(O, i + 1, PB, PA, s[i])
        end
    end
    O
end

# true if v is a permutation of 1:D (D = length(v))
function _is_permutation_1toD!(seen, v)
    eltype(v) <: Integer || return false
    fill!(seen, false)
    D = length(v)
    @inbounds for x in v
        (1 <= x <= D && !seen[x]) || return false
        seen[x] = true
    end
    true
end

# row i of O keeps P1[1:k] and is completed with the genes of P2 not in P1[1:k],
# in the order they appear in P2.
function _order_fill!(O, i, P1, P2, k, seen)
    fill!(seen, false)
    @inbounds for j in 1:k
        seen[P1[j]] = true
    end
    c = k
    @inbounds for x in P2
        seen[x] && continue
        c += 1
        O[i, c] = x
    end
    O
end

function _order_fill_generic!(O, i, P1, P2, k)
    tail = setdiff(P2, view(P1, 1:k))
    if length(tail) != size(O, 2) - k
        throw(ArgumentError(
            "OrderCrossover requires each row to be a permutation of the same values.",
        ))
    end
    O[i, k + 1:end] = tail
    O
end
