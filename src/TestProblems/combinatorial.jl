import ..PermutationSpace,  ..BitArraySpace

"""
    knapsack(profit, weight, capacity; encoding=:permutation)

Return the Knapsack problem regarding the provided parameters. 
"""
function knapsack(profit, weight, capacity; encoding=:permutation)
    function _f(items)
        p = w = 0.0
        for item in items
            w += weight[item]
            if w > capacity
                return -p
            end
            p += profit[item]
        end
        -p
    end
    if encoding === :permutation
        f = _f
        search_space = PermutationSpace(length(profit))
    elseif encoding === :binary
        f = x -> _f(findall(x))
        search_space = BitArraySpace(length(profit))
    else
        error("Encoding $encoding is not valid.")
    end
    # TODO Try to handle optimums when provided when provided
    optimum = nothing

    return f, search_space, optimum
end
