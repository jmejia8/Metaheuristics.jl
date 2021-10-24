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

function DTLZ_hyperplane!(fx, x)
    m = length(fx)
    for i in 1:m
        fx[i] *= prod(x[1:m-i])
        if i > 1
            fx[i] *= 1 - x[1+m - i]
        end
    end
    fx
end

function DTLZ_g1(x, m)
    y = view(x, m:length(x)) .- 0.5
    return 100*( length(y) + sum( y.^2 - cos.(20π*y) ))
end

function DTLZ_g2(x, m)
    return sum( (view(x, m:length(x)) .- 0.5).^2 )
end


function DTLZ1_f(x, m = 3)
    g = DTLZ_g1(x, m)
    D = length(x)

    fx = fill(0.5*(1 + g), m)
    DTLZ_hyperplane!(fx, x)

    return fx, [0.0], [0.0]
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

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = 0.5ref_dirs
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return DTLZ1_f, bounds, pareto_set
end


function DTLZ2_f(x, m = 3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x)

    return fx, [0.0], [0.0]
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

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return DTLZ2_f, bounds, pareto_set
end


function DTLZ3_f(x,m=3)
    g = DTLZ_g1(x, m)
    fx = fill(1 + g, m)
    DTLZ_hypersphere!(fx, x)

    return fx, [0.0], [0.0]
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

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return DTLZ3_f, bounds, pareto_set
end


function DTLZ4_f(x,m=3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x; α = 100)

    return fx, [0.0], [0.0]
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

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return DTLZ4_f, bounds, pareto_set
end

function DTLZ5_f(x,m = 3)
    g = DTLZ_g2(x, m)
    fθ = fill(1.0 + g, m)
    θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:m-1] )
    θ[1] = x[1]


    DTLZ_hypersphere!(fθ, θ)

    return fθ, [0.0], [0.0]
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

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    if m == 3
        ref_dirs = gen_ref_dirs(2, 100)
        X = fill(0.5, length(ref_dirs), D)

        for i in eachindex(ref_dirs)
            X[i,1:2] = ref_dirs[i]
        end
        pareto_set = [ generateChild(X[i,:], DTLZ5_f(X[i,:])) for i in eachindex(ref_dirs) ]

    else
        pareto_set = []
    end

    return DTLZ5_f, bounds, pareto_set
end


function DTLZ6_f(x,m=3)
    g = sum(x[m:end] .^ 0.1)
    fθ = fill(1.0 + g, m)
    θ = @. 1 / (2*(1 + g)) * (1 + 2g * x[1:m-1] )
    θ[1] = x[1]


    DTLZ_hypersphere!(fθ, θ)

    return fθ, [0.0], [0.0]
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

    D = 10 

    bounds = Array([zeros(D) ones(D)]')

    if m == 3
        ref_dirs = gen_ref_dirs(2, 100)
        X = fill(0.0, length(ref_dirs), D)

        for i in eachindex(ref_dirs)
            X[i,1:2] = ref_dirs[i]
        end
        pareto_set = [ generateChild(X[i,:], DTLZ6_f(X[i,:])) for i in eachindex(ref_dirs) ]

    else
        pareto_set = []
    end

    return DTLZ6_f, bounds, pareto_set
end

#=
####################################################################################
####################################################################################
####################################################################################
#                Constrained DTLZ
####################################################################################
####################################################################################
####################################################################################
=#

function C1_DTLZ1_f(x, m = 3)
    g = DTLZ_g1(x, m)
    D = length(x)

    fx = fill(0.5*(1 + g), m)
    DTLZ_hyperplane!(fx, x)


    c = fx[end]/0.6 + sum(fx[1:end-1]/0.5) - 1;
    return fx, [c], [0.0]
end


"""
    C1_DTLZ1(m = 3, ref_dirs = gen_ref_dirs(m, 12))

C1_DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- convex
- multifrontal
- constraints type 1
"""
function C1_DTLZ1(m=3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = 0.5ref_dirs
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C1_DTLZ1_f, bounds, pareto_set
end


function C1_DTLZ3_f(x,m=3)
    g = DTLZ_g1(x, m)
    fx = fill(1 + g, m)
    DTLZ_hypersphere!(fx, x)

    if m == 2
        r = 6
    elseif m <= 3
        r = 9
    elseif m <= 8
        r = 12.5
    else
        r = 15
    end
    c = -(sum(fx .^ 2) - 16) * (sum(fx .^ 2) - r^2);

    return fx, [c], [0.0]
end

"""
    C1_DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ3 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- multifrontal
- constraints type 1
"""
function C1_DTLZ3(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C1_DTLZ3_f, bounds, pareto_set
end


function C2_DTLZ2_f(x, m = 3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x)

    if m == 3
        r = 0.4
    else
        r = 0.5
    end
    

    p1 = maximum( i -> (fx[i] .- 1).^2 + sum(fx[1:m].^2 .- fx[i]^2) - r^2,1:m)
    p2 = sum((fx .- 1 / sqrt(m)).^2) - r^2
    return fx, [-max(p1, p2)], [0.0]
end

"""
    C2_DTLZ2(m = 3, ref_dirs = gen_ref_dirs(m, 12))

DTLZ2 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- unifrontal
- contraints type 2
"""
function C2_DTLZ2(m=3, ref_dirs = gen_ref_dirs(m, 12))

    D = 10 + m - 1

    bounds = Array([zeros(D) ones(D)]')

    pf = generic_sphere(ref_dirs)
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C2_DTLZ2_f, bounds, pareto_set
end



function C3_DTLZ1_f(x, m = 3)
    g = DTLZ_g1(x, m)
    D = length(x)

    fx = fill(0.5*(1 + g), m)
    DTLZ_hyperplane!(fx, x)


    c = [sum(fx[j] .+ fx/0.5 .- fx[j]/0.5) - 1 for j in 1:m]
    return fx, -c, [0.0]
end

#=
"""
    C3_DTLZ1(m = 3, ref_dirs = gen_ref_dirs(m, 12))

C3_DTLZ1 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- convex
- multifrontal
- constraints type 3
"""
function C3_DTLZ1(m=3, ref_dirs = gen_ref_dirs(m, 12))

    D = m + 4

    bounds = Array([zeros(D) ones(D)]')

    @warn "C3_DTLZ1 is under development"

    pf = ref_dirs./([(0.5r -6*maximum(r)) for r in ref_dirs])
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C3_DTLZ1_f, bounds, pareto_set
end
=#




function C3_DTLZ4_f(x,m=3)
    g = DTLZ_g2(x, m)
    fx = fill(1.0 + g, m)
    DTLZ_hypersphere!(fx, x; α = 100)

    c = [fx[j]^2/4  + sum(fx.^2 .- fx[j]^2) - 1 for j in 1:m]
    return fx, -c, [0.0]
end

"""
    C3_DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))

C3_DTLZ4 returns `(f::function, bounds::Matrix{Float64}, pareto_set::Array{xFgh_indiv})`
where `f` is the objective function and `pareto_set` is an array with optimal Pareto solutions
with `n_solutions`.

### Parameters
- `m` number of objective functions
- `ref_dirs` number of Pareto solutions (default: Das and Dennis' method).

Main properties:
- nonconvex
- unifrontal
- constraints type 3
"""
function C3_DTLZ4(m = 3, ref_dirs = gen_ref_dirs(m, 12))

    D = m + 4

    bounds = Array([zeros(D) ones(D)]')

    pf = ref_dirs./([sqrt.(sum(r.^2)-3/4*maximum(r.^2)) for r in ref_dirs])
    
    pareto_set = [ generateChild(zeros(0), (fx, [0.0], [0.0])) for fx in pf ]

    return C3_DTLZ4_f, bounds, pareto_set
end
