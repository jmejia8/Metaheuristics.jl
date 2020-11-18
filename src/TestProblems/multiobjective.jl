"""
convex
"""
function ZDT1(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return ( [x[1], gx*(1 - sqrt(x[1] / gx))  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    X = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in X ]
     
    return f, bounds, pareto_set
end

"""
nonconvex
"""
function ZDT2(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return ( [x[1], gx*(1 - (x[1] / gx)^2)  ], [0.0], [0.0] )
    end
    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(d - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]
     
    return f, bounds, pareto_set
end

"""
convex disconected
"""
function ZDT3(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        a = x[1] / gx
        return ( [x[1], gx*(1 - sqrt(a) - a*sin(10π*x[1]) )  ], [0.0], [0.0] )
    end
    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]
     
    return f, bounds, pareto_set

end

"""
nonconvex
"""
function ZDT4(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 10*(length(x)-1) + sum( x[2:end].^2 - 10cos.(4π*x[2:end])) 
        return ( [x[1], gx*(1 - (x[1] / gx)^2)  ], [0.0], [0.0] )
    end
    bounds = Array([-5zeros(D) 5ones(D)]')
    bounds[:,1] = [1, 1.0]

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]
     
    return f, bounds, pareto_set
end

"""
nonconvex nonunifromly spaced
"""
function ZDT6(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )^(0.25)
        ff1 = 1 - exp(-4x[1])*sin(6π*x[1])^6 
        return ( [ ff1 , gx*(1 - (ff1 / gx)^2)  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]
     
    return f, bounds, pareto_set
end

