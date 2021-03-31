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

function constrained3()
    f(y) = begin
        x = zeros(2)
        return x[1] + 2x[2] + y[1] + y[2] + 2y[3],
                    [0.0], # g
                    [ 
                     - y[1] + y[2] + y[3] + y[4] - 1,
                     2x[1] - y[1] + 2y[2] - 0.5y[3] + y[5] - 1,
                     2x[2] + 2y[1] - y[2] - 0.5y[3] + y[6] - 1
                    ]
    end


    bounds = Array([zeros(6) 100.0ones(6)]')
    y = zeros(6)
    return f, bounds, [generateChild(y, f(y)) ]
end
