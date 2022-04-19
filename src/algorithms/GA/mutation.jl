"""
    BitFlipMutation(;p = 1e-2)

Flip each bit with probability `p`.
"""
struct BitFlipMutation
    p::Float64
    BitFlipMutation(;p = 1e-2) = new(p)
end

function mutation!(Q, parameters::BitFlipMutation)
    mask = rand(size(Q)...) .< parameters.p
    Q[mask] = .!Q[mask]
    Q
end

"""
    SlightMutation

> Fogel, D. B. (1988). An evolutionary approach to the traveling salesman problem.
> Biological Cybernetics, 60(2), 139-144.
"""
struct SlightMutation end

function mutation!(Q, parameters::SlightMutation)
    N, D = size(Q)
    k = rand(1:D, N)
    s = rand(1:D, N)
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
