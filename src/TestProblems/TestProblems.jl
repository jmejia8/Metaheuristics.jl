module TestProblems

import ..generateChild, ..gen_ref_dirs, ..norm, ..get_non_dominated_solutions

include("box-constrained.jl")
include("constrained.jl")
include("multiobjective.jl")
include("combinatorial.jl")

const _problems_dict = Dict(
                      :sphere => sphere,
                      :discus => discus,
                      :rastrigin => rastrigin,
                      #################################
                      :ZDT1 => ZDT1,
                      :ZDT2 => ZDT2,
                      :ZDT3 => ZDT3,
                      :ZDT4 => ZDT4,
                      :ZDT6 => ZDT6,
                      #################################
                      :constrained1 => constrained1,
                      :constrained2 => constrained2,
                      :constrained3 => constrained3,
                      #################################
                      :MTP => MTP,
                      #################################
                      :DTLZ1 => DTLZ1,
                      :DTLZ2 => DTLZ2,
                      :DTLZ3 => DTLZ3,
                      :DTLZ4 => DTLZ4,
                      :DTLZ5 => DTLZ5,
                      :DTLZ6 => DTLZ6,
                      #################################
                      :C1_DTLZ1 => C1_DTLZ1,
                      :C1_DTLZ3 => C1_DTLZ3,
                      :C2_DTLZ2 => C2_DTLZ2,
                      :C3_DTLZ4 => C3_DTLZ4,
                      #################################
                      :knapsack => knapsack,
                     )

"""
    get_problem(problem)

Returns a 3-tuple with the objective function, the bounds and 100 Pareto solutions for
multi-objective optimization problems or the optimal solutions for (box)constrained optimization
problems.

Here, `problem` can be one of the following symbols:

Single-objective:
- `:sphere`
- `:discus`
- `:rastrigin`

Multi-objective:
- `:ZDT1`
- `:ZDT2`
- `:ZDT3`
- `:ZDT4`
- `:ZDT6`


Many-objective:
- `:DTLZ1`
- `:DTLZ2`
- `:DTLZ3`
- `:DTLZ4`
- `:DTLZ5`
- `:DTLZ6`

## Example

```jldoctest
julia> import Metaheuristics: TestProblems, optimize

julia> f, bounds, pareto_solutions = TestProblems.get_problem(:ZDT3);


julia> bounds
2×30 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  …  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0
 1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0     1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0

julia> pareto_solutions
                           F space
          ┌────────────────────────────────────────┐ 
        1 │⢅⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠈⢢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠙⠒⠀⠀⠀⠀⡆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠘⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠱⠄⠀⠀⠀⠀⠀⠀⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢱⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   f_2    │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠬⡦⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠘⠂⠀⠀⠀⠀⠀⠀⢢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠸⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢧⡀⠀⠀⠀⠀⠀⠀⢀⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⢸⠀⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⡆⠀⠀│ 
          │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠙⠀⠀│ 
       -1 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
          └────────────────────────────────────────┘ 
          0                                      0.9
                             f_1

```
"""
function get_problem(problem)
    # TODO improve this part
    if problem in keys(_problems_dict)
        return _problems_dict[problem]()
    end

    _p = keys(_problems_dict)
    error("Problem $problem is not defined. Implemented problems:\n$_p")
end

end
