"""
    SlightMutation

> Fogel, D. B. (1988). An evolutionary approach to the traveling salesman problem.
> Biological Cybernetics, 60(2), 139-144.
"""
struct SlightMutation
    rng
end

SlightMutation(;rng = default_rng_mh()) = SlightMutation(rng)

function mutation!(Q, parameters::SlightMutation)
    N, D = size(Q)
    k = rand(parameters.rng, 1:D, N)
    s = rand(parameters.rng, 1:D, N)
    for i in 1:N
        s[i] == k[i] && continue
        if s[i] < k[i]
            c = vcat(1:s[i]-1, k[i], s[i]:k[i]-1, k[i]+1:D)
        else s[i] > k[i]
            c = vcat(1:k[i]-1, k[i]+1:s[i]-1, k[i], s[i]:D)
        end
        Q[i,:] = Q[i, c]
    end
    Q
end

