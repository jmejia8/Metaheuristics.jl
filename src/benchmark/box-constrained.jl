function sphere(D=10)
    # Objective function
    f(x) = sum(x.*x)

    bounds = Array([-100.0ones(D) 100.0ones(D)]')

	return f, bounds

end

function discus(D = 10)
    # Objective function
    f(x) = 1e6x[1].^2 + sum(x[2:end] .^2)

    bounds = Array([-10.0ones(D) 10.0ones(D)]')
	
    return f, bounds
end

function rastrigin(D = 2)
    
    # Objective function
    f(x) = 10D+ sum(x.*x - 10cos.(2π*x))

    bounds = Array([-5.0ones(D) 5.0ones(D)]')

    return f, bounds

end

