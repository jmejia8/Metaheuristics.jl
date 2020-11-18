"""
    sphere(D)

The well-known D-dimensional Sphere function.
"""
function sphere(D=10)
    # Objective function
    f(x) = sum(x.*x)

    bounds = Array([-100.0ones(D) 100.0ones(D)]')

    x = zeros(D)
    return f, bounds, [generateChild(x, f(x)) ]

end

"""
    discus(D)

The well-known D-dimensional Discus function.
"""
function discus(D = 10)
    # Objective function
    f(x) = 1e6x[1].^2 + sum(x[2:end] .^2)

    bounds = Array([-10.0ones(D) 10.0ones(D)]')
	
    x = zeros(D)
    return f, bounds, [generateChild(x, f(x))] 
end


"""
    rastrigin(D)

The well-known D-dimensional Rastrigin function.
"""
function rastrigin(D = 10)
    
    # Objective function
    f(x) = 10D+ sum(x.*x - 10cos.(2π*x))

    bounds = Array([-5.0ones(D) 5.0ones(D)]')

    x = zeros(D)
    return f, bounds, [generateChild(x, f(x)) ]

end

