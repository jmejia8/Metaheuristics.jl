"""
convex
"""
function ZDT1()
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return ( [x[1], gx*(1 - sqrt(x[1] / gx))  ], [0.0], [0.0] )
    end
    bounds = Array([zeros(30) ones(30)]')
    return f, bounds
end

"""
nonconvex
"""
function ZDT2()
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return ( [x[1], gx*(1 - (x[1] / gx)^2)  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(30) ones(30)]')
    return f, bounds
end

"""
convex disconected
"""
function ZDT3()
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        a = x[1] / gx
        return ( [x[1], gx*(1 - sqrt(a) - a*sin(10π*x[1]) )  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(30) ones(30)]')
    return f, bounds

end

"""
nonconvex
"""
function ZDT4()
    f(x) = begin
        gx = 1.0 + 10*(length(x)-1) + sum( x[2:end].^2 - 10cos.(4π*x[2:end])) 
        return ( [x[1], gx*(1 - (x[1] / gx)^2)  ], [0.0], [0.0] )
    end
    bounds = Array([-5zeros(10) 5ones(10)]')
    bounds[:,1] = [1, 1.0]
    return f, bounds
end

"""
nonconvex nonunifromly spaced
"""
function ZDT6()
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )^(0.25)
        ff1 = 1 - exp(-4x[1])*sin(6π*x[1])^6 
        return ( [ ff1 , gx*(1 - (ff1 / gx)^2)  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(10) ones(10)]')
    return f, bounds
end

