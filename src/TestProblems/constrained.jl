function constrained1(D = 5)
    # Objective function
    f(x) = (sum((x .- 1).^2), [sum((x .- 1).^2) - 4, sum(sin.(x .- 1)) - 1], [(x[1] - x[2])^2])


    bounds = Array([-10.0ones(D) 10.0ones(D)]')

    x = zeros(D)
    return f, bounds, [generateChild(x, f(x)) ]
end


function constrained2(D = 5)

    f(x) = (sum((x .- 1).^2), [sum((x .- 1).^2) - 4, sum(sin.(x .- 1)) - 1, (x[1] - x[2])^2 - 6])


    bounds = Array([-10.0ones(D) 10.0ones(D)]')

    x = zeros(D)
    return f, bounds, [generateChild(x, f(x)) ]
    
end

