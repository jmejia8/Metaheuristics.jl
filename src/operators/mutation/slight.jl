"""
    SlightMutation

> Fogel, D. B. (1988). An evolutionary approach to the traveling salesman problem.
> Biological Cybernetics, 60(2), 139-144.
"""
struct SlightMutation{R}
    rng::R
end

SlightMutation(; rng = default_rng_mh()) = SlightMutation(rng)

function mutation!(Q, parameters::SlightMutation)
    N, D = size(Q)
    k = rand(parameters.rng, 1:D, N)
    s = rand(parameters.rng, 1:D, N)
    # move the gene at k[i] to position s[i] (rotate the segment between them by one)
    @inbounds for i in 1:N
        ki, si = k[i], s[i]
        ki == si && continue
        t = Q[i, ki]
        if si < ki
            for j in ki:-1:si + 1
                Q[i, j] = Q[i, j - 1]
            end
            Q[i, si] = t
        else
            for j in ki:si - 2
                Q[i, j] = Q[i, j + 1]
            end
            Q[i, si - 1] = t
        end
    end
    Q
end
