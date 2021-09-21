generic_sphere(ref_dirs) = ref_dirs ./ ([norm(v) for v in ref_dirs])

function DTLZ_hypersphere!(fx, x; α = 1)
    m = length(fx)
    for i = 1:m
        if i < m
            fx[i] *= prod(cos.((π*0.5) * x[1:m-i] .^ α ))
        end
        if i > 1
            fx[i] *= sin((π*0.5) * x[1+m - i] .^ α)
        end
    end

    return fx
end

function DTLZ_g1(x, m)
    y = view(x, m:length(x)) .- 0.5
    return 100*( length(y) + sum( y.^2 - cos.(20π*y) ))
end

function DTLZ_g2(x, m)
    return sum( (view(x, m:length(x)) .- 0.5).^2 )
end



"""
    DTLZ1(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- convex
- multifrontal
"""
function DTLZ1(m=3, ref_dirs = gen_ref_dirs(m, 12))
    f(x,m=m) = begin
        g = DTLZ_g1(x, m)
        D = length(x)

        fx = fill(0.5*(1 + g), m)
        for i in 1:m
            fx[i] *= prod(x[1:m-i])
            if i > 1
                fx[i] *= 1 - x[1+m - i]
            end
        end

        return fx, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = 0.5ref_dirs
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return f, bounds, pareto_set
end


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
        g = DTLZ_g2(x, m)
        fx = fill(1.0 + g, m)
        DTLZ_hypersphere!(fx, x)

        return fx, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return f, bounds, pareto_set
end


"""
    DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ3 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- multifrontal
"""
function DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))
    f(x,m=m) = begin
        g = DTLZ_g1(x, m)
        fx = fill(1 + g, m)
        DTLZ_hypersphere!(fx, x)

        return fx, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return f, bounds, pareto_set
end


"""
    DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ4 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- unifrontal
"""
function DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))
    f(x,m=m) = begin
        g = DTLZ_g2(x, m)
        fx = fill(1.0 + g, m)
        DTLZ_hypersphere!(fx, x; α = 100)

        return fx, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return f, bounds, pareto_set
end



"""
    DTLZ5(m = 3)

DTLZ5 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions

"""
function DTLZ5(m = 3)
    f(x,m=m) = begin
        g = DTLZ_g2(x, m)
        fθ = fill(1.0 + g, m)
        θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:m-1] )
        θ[1] = x[1]


        DTLZ_hypersphere!(fθ, θ)

        return fθ, [0.0], [0.0]
    end

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    if m == 3
        ref_dirs = gen_ref_dirs(2, 100)
        X = fill(0.5, length(ref_dirs), D)

        for i in eachindex(ref_dirs)
            X[i,1:2] = ref_dirs[i]
        end
        pareto_set = [ generateChild(X[i,:], f(X[i,:])) for i in eachindex(ref_dirs) ]
        
    else
        pareto_set = []
    end
    
    return f, bounds, pareto_set
end


"""
    DTLZ6(m = 3)

DTLZ6 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions

"""
function DTLZ6(m = 3)
    f(x,m=m) = begin
        g = sum(x[m:end] .^ 0.1)
        fθ = fill(1.0 + g, m)
        θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:m-1] )
        θ[1] = x[1]


        DTLZ_hypersphere!(fθ, θ)

        return fθ, [0.0], [0.0]
    end

    D = 10 

    bounds = Array([zeros(D) ones(D)]')

    if m == 3
        ref_dirs = gen_ref_dirs(2, 100)
        X = fill(0.0, length(ref_dirs), D)

        for i in eachindex(ref_dirs)
            X[i,1:2] = ref_dirs[i]
        end
        pareto_set = [ generateChild(X[i,:], f(X[i,:])) for i in eachindex(ref_dirs) ]
        
    else
        pareto_set = []
    end

    return f, bounds, pareto_set
end
