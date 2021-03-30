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

# for Deb rules
function countViolations(g::Vector, h::Vector)
    sum_g = 0
    sum_h = 0

    for i = 1:length(g)
        if g[i] > 0
        sum_g += 1  end
    end

    for i = 1:length(h)
        if h[i] != 0.0
        sum_h += 1  end
    end

    return sum_g + sum_h
end

function isfeasible(element::xf_indiv)
    return true
end


function isfeasible(element::xfgh_indiv)
    return countViolations(element.g, element.h) == 0
end

