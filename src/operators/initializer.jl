abstract type AbstractInitializer end

struct RandomInBounds <: AbstractInitializer
    N::Int
end

"""
    RandomInBounds

Initialize `N` solutions with random values in bounds. Suitable for
integer and real coded problems.
"""
RandomInBounds(;N = 0) = RandomInBounds(N)

"""
    RandomBinary(;N)

Create random binary individuals.
"""
struct RandomBinary <: AbstractInitializer
    N::Int
    RandomInBounds(;N = 0) = new(N)
end




"""
    RandomPermutation(;N)

Create individuals in random permutations.
"""
struct RandomPermutation <: AbstractInitializer
    N::Int
    RandomPermutation(;N = 0) = new(N)
end


"""
    LatinHypercubeSampling(nsamples, dim; iterations)

Create `N` solutions within a Latin Hypercube sample in bounds with dim.

### Example 

```julia-repl
julia> sample(LatinHypercubeSampling(10,2))
10×2 Matrix{Float64}:
 0.0705631  0.795046
 0.7127     0.0443734
 0.118018   0.114347
 0.48839    0.903396
 0.342403   0.470998
 0.606461   0.275709
 0.880482   0.89515
 0.206142   0.321041
 0.963978   0.527518
 0.525742   0.600209

julia> sample(LatinHypercubeSampling(10,2), [-10 -10;10 10.0])
10×2 Matrix{Float64}:
 -7.81644   -2.34461
  0.505902   0.749366
  3.90738   -8.57816
 -2.05837    9.803
  5.62434    6.82463
 -9.34437    2.72363
  6.43987   -1.74596
 -1.3162    -4.50273
  9.45114   -7.13632
 -4.71696    5.0381
```
"""
struct LatinHypercubeSampling <: AbstractInitializer
    N::Int
    dim::Int
    iterations::Int
    LatinHypercubeSampling(nsamples, dim;iterations=25) = new(nsamples,dim,iterations)
end

"""
    Grid(npartitions, dim)

Parameters to generate a grid with `npartitions` in a space with `dim` dimensions.

### Example 

```julia-repl
julia> sample(Grid(5,2))
25×2 Matrix{Float64}:
 0.0   0.0
 0.25  0.0
 0.5   0.0
 0.75  0.0
 ⋮     
 0.5   1.0
 0.75  1.0
 1.0   1.0

julia> sample(Grid(5,2), [-1 -1; 1 1.])
25×2 Matrix{Float64}:
 -1.0  -1.0
 -0.5  -1.0
  0.0  -1.0
  0.5  -1.0
  ⋮    
  0.0   1.0
  0.5   1.0
  1.0   1.0
```

Note that the sample is with size `npartitions^(dim)`.
"""
struct Grid <: AbstractInitializer
    npartitions::Int
    dim::Int
end

function gen_initial_state(problem,parameters::RandomPermutation,information,options)

    D = getdim(problem)
    X = zeros(Int, parameters.N, D)
    N = parameters.N
    for i in 1:parameters.N
        X[i,:] = shuffle(1:D)
    end

    if problem.parallel_evaluation
        population = create_solutions(X, problem; ε=options.h_tol)
    else
        population = [ create_solution(X[i,:], problem; ε=options.h_tol) for i in 1:N]
    end 

    best_solution = get_best(population)

    State(best_solution, population; f_calls = length(population), iteration=1)
end

function initialize!(
        status,
        parameters::AbstractInitializer,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    gen_initial_state(problem,parameters,information,options, status)
end

include("sample.jl")

