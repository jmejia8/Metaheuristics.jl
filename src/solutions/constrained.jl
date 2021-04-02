# COP functions
function violationsSum(g::Vector, h::Vector; ε =0.0)
    sum_g = sum(max.(0.0, g))
    sum_h = 0.0


    for i = 1:length(h)
        if !isapprox(h[i], 0.0, atol=ε)
            sum_h += abs(h[i])
        end
    end

    return (sum_g/length(g)) + (sum_h / length(h))
end

