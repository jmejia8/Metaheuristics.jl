"""
    ZDT1(D, n_solutions)

ZDT1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `D` number of variables (dimension)
- `n_solutions` number of pareto solutions.

Main properties:
- convex
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
    ZDT2(D, n_solutions)

ZDT2 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `D` number of variables (dimension)
- `n_solutions` number of pareto solutions.

Main properties:
- nonconvex
"""
function ZDT2(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        return ( [x[1], gx*(1 - (x[1] / gx)^2)  ], [0.0], [0.0] )
    end
    bounds = Array([zeros(D) ones(D)]')

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]

    return f, bounds, pareto_set
end

"""
    ZDT3(D, n_solutions)

ZDT3 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `D` number of variables (dimension)
- `n_solutions` number of pareto solutions.

Main properties:
- convex disconected
"""
function ZDT3(D = 30, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1) )
        a = x[1] / gx
        return ( [x[1], gx*(1 - sqrt(a) - a*sin(10π*x[1]) )  ], [0.0], [0.0] )
    end
    bounds = Array([zeros(D) ones(D)]')

    regions = [ 0 0.0830015349
                0.182228780 0.2577623634
                0.4093136748 0.4538821041
                0.6183967944 0.6525117038
                0.8233317983 0.8518328654]

    n = n_solutions ÷ size(regions, 1)
    x = Float64[]
    for i in 1:size(regions, 1)
        x = vcat(x, range(regions[i,1], regions[i,2], length=n))
    end

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]

    return f, bounds, get_non_dominated_solutions(pareto_set)

end

"""
    ZDT4(D, n_solutions)

ZDT4 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `D` number of variables (dimension)
- `n_solutions` number of pareto solutions.

Main properties:
- nonconvex
"""
function ZDT4(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 10*(length(x)-1) + sum( x[2:end].^2 - 10cos.(4π*x[2:end]))
        return ( [x[1], gx*(1 - sqrt(x[1] / gx))  ], [0.0], [0.0] )
    end
    bounds = Array([-5zeros(D) 5ones(D)]')
    bounds[:,1] = [0, 1.0]

    x = range(0, 1, length=n_solutions)

    x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    pareto_set = [ generateChild(x, f(x)) for x in x ]

    return f, bounds, pareto_set
end

"""
    ZDT6(D, n_solutions)

ZDT6 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `D` number of variables (dimension)
- `n_solutions` number of Pareto solutions.

Main properties:
- nonconvex
- non-uniformly spaced
"""
function ZDT6(D = 10, n_solutions = 100)
    f(x) = begin
        gx = 1.0 + 9.0 * ( sum(x[2:end]) / (length(x)-1.0) )^(0.25)
        ff1 = 1.0 - exp(-4.0x[1])*sin(6.0π*x[1])^6
        return ( [ ff1 , gx*(1.0 - (ff1 / gx)^2)  ], [0.0], [0.0] )
    end

    bounds = Array([zeros(D) ones(D)]')

    #x = range(0, 1, length=n_solutions)

    #x = [vcat(x[i], zeros(D - 1)) for i in 1:n_solutions]
    xx = range(0.2807753191, 1, length=100)
    yy = 1 .- (xx).^2
    pareto_set = [ generateChild(zeros(0), ([xx[i], yy[i]], [0.0], [0.0])) for i in 1:length(xx) ]

    return f, bounds, pareto_set
end

function MTP(D = 10, n_solutions = 100)
    f(x) = begin
        Q = 1.0 .+ sum((x[3:end] .- 0.8).^2)
        fx = [x[1], x[2]] .* Q
        g = [ (x[1] - 1)^2 + (x[2] - 1)^2 - 1 ]
        h = [0.0]
        return fx, g, h
    end

    bounds = Array([ zeros(D) 2.0ones(D)]')
    θ = range(π, 1.5π, length=n_solutions)

    x = 1.0 .+ cos.(θ)
    y = 1.0 .+ sin.(θ)

    pareto_set = [ generateChild(zeros(0), ([x[i], y[i]], [0.0], [0.0])) for i in 1:length(θ) ]

    return f, bounds, pareto_set


end

generic_sphere(ref_dirs) = ref_dirs ./ ([norm(v) for v in ref_dirs])

"""
    DTLZ2(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ2 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- unifrontal
"""
function DTLZ2(m=3, ref_dirs = gen_ref_dirs(m, 12))
    f(x,m=m) = begin
        g = sum( (view(x, m:length(x)) .- 0.5).^2 )
        fx = zeros(m)
        for i = 1:m
            fx[i] = (1.0 + g)
            if i < m
                fx[i] *= prod(cos.(x[1:m-i] * (π*0.5) ))
            end
            if i > 1
                fx[i] *= sin( x[1+m - i] * (π*0.5) )
            end
        end
        return fx, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return f, bounds, pareto_set
end
