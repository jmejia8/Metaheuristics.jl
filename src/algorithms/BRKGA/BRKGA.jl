"""
    BRKGA(num_elites = 20, num_mutants = 10, num_offsprings = 70, bias = 0.7)

Biased Random Key Genetic Algorithm (BRKGA).

### Example

```julia-repl
julia> target_perm = collect(reverse(1:10))
10-element Vector{Int64}:
 10
  9
  8
  7
  6
  5
  4
  3
  2
  1

julia> decode(rk) = sortperm(rk);

julia> f(rk) = sum(abs.(decode(rk) - target_perm));

julia> res = optimize(f, [zeros(10) ones(10)], BRKGA(num_elites=70))
Optimization Result
===================
  Iteration:       22
  Minimum:         0
  Minimizer:       [0.938595, 0.851247, 0.823736, …, 0.375321]
  Function calls:  2200
  Total time:      0.0238 s
  Stop reason:     Due to Convergence Termination criterion.

julia> decode(minimizer(res))
10-element Vector{Int64}:
 10
  9
  8
  7
  6
  5
  4
  3
  2
  1
```
"""
function BRKGA(;
        num_elites = 20,
        num_mutants = 10,
        num_offsprings = 70,
        N = num_elites + num_mutants + num_offsprings,
        bias = 0.7,
        kargs...
    )

    initializer = RandomInBounds(;N)
    selection   = BiasedSelection(num_elites, num_offsprings)
    crossover   = BinomialCrossover(p = bias, n_offsprings = 1)
    mutation    = InsertRandomMutation(num_mutants)
    environmental_selection = ElitistReplacement()

    GA(;initializer,
       selection,
       crossover,
       mutation,
       environmental_selection,
       kargs...
      )
end

