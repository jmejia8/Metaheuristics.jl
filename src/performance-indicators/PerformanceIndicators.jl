
"""
	PerformanceIndicators

This module includes performance indicators to assess evolutionary multi-objective
optimization algorithms.

- `gd` Generational Distance
- `igd` Inverted Generational Distance
- `gd_plus` Generational Distance plus
- `igd_plus` Inverted Generational Distance plus

### Example

```jldoctest
julia> import Metaheuristics: PerformanceIndicators, TestProblems

julia> A = [ collect(1:3) collect(1:3) ]
3×2 Array{Int64,2}:
 1  1
 2  2
 3  3

julia> B = A .- 1
3×2 Array{Int64,2}:
 0  0
 1  1
 2  2

julia> PerformanceIndicators.gd(A, B)
0.47140452079103173

julia> f, bounds, front = TestProblems.get_problem(:ZDT1);

julia> front
                          F space
         ┌────────────────────────────────────────┐ 
       1 │⠁⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠈⠄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠈⢆⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠈⠢⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠉⠢⡄⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
   f_2   │⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠲⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠒⢤⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠙⠢⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠑⠢⢄⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠉⠢⠤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│ 
         │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠑⠢⢤⣀⠀⠀⠀⠀⠀│ 
       0 │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠉⠒⠢⢄⣀│ 
         └────────────────────────────────────────┘ 
         0                                        1
                            f_1

julia> PerformanceIndicators.igd_plus(front, front)
0.0
```


"""
module PerformanceIndicators

import ..State, ..fvals, ..norm, ..xFgh_indiv, ..fval, ..compare, ..mean
import ..get_non_dominated_solutions

"""
	generational_distance(front, true_pareto_front; p = 1, inverted = false, plus=false)

Returns the Generational Distance.

### Parameters

- `front` and `true_pareto_front` can be:
	- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
	- `State`
	- `Array{xFgh_indiv}` (usually `State.population`)
- `true_pareto_front` is a `M×m` matrix.
- `p` is the power in for the ‖⋅‖_p
- `plus` if true then computes the GD+
- `inverted` if true then computes IGD
"""
function generational_distance(front, true_pareto_front; p = 1, inverted = false, plus=false)

	if inverted
		tmp = true_pareto_front
		true_pareto_front = front
		front = tmp
	end

	distances = zeros(size(true_pareto_front, 1))

	s = 0.0
	for i in 1:size(front,1)	
		x = view(front, i, :)
		for j = 1:size(true_pareto_front, 1)
			y = view(true_pareto_front, j, :)
			if plus	&& inverted
				distances[j] = norm( max.(y - x, 0) )		
			elseif plus && !inverted
				distances[j] = norm( max.(x - y, 0) )		
			else
				distances[j] = norm( y - x)		
			end
			
		end	

		s += minimum(distances) ^ p
	end

	return (s^(1/p)) /  size(front, 1)
	
end


function generational_distance(front::Union{Array{xFgh_indiv}, State},
						       true_pareto_front::Union{Array{xFgh_indiv}, State};
								p = 1,
								inverted = false,
								plus = false
		)
	generational_distance(fvals(front), fvals(true_pareto_front); p=p,inverted=inverted,plus=plus)
end


function generational_distance(front::Union{Array{xFgh_indiv}, State},
						       true_pareto_front;
								p = 1,
								inverted = false,
								plus = false
		)
	generational_distance(fvals(front), (true_pareto_front); p=p,inverted=inverted,plus=plus)
end


"""
	gd(front, true_pareto_front; p = 1)

Returns the Generational Distance.

### Parameters

`front` and `true_pareto_front` can be:
	- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
	- `State`
	- `Array{xFgh_indiv}` (usually `State.population`)
"""
gd(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, p=p)


"""
	igd(front, true_pareto_front; p = 1)

Returns the Inverted Generational Distance.

### Parameters

`front` and `true_pareto_front` can be:
	- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
	- `State`
	- `Array{xFgh_indiv}` (usually `State.population`)
"""
igd(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, inverted=true, p=p)

"""
	gd_plus(front, true_pareto_front; p = 1)

Returns the Generational Distance Plus.

### Parameters

`front` and `true_pareto_front` can be:
	- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
	- `State`
	- `Array{xFgh_indiv}` (usually `State.population`)
"""
gd_plus(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, p=p, plus=true)

"""
	igd_plus(front, true_pareto_front; p = 1)

Returns the Inverted Generational Distance Plus.

### Parameters

`front` and `true_pareto_front` can be:
	- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
	- `State`
	- `Array{xFgh_indiv}` (usually `State.population`)
"""
igd_plus(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, inverted=true, p=p, plus=true)


"""
	spacing(A)
Computes the Schott spacing indicator. `spacing(A) == 0` means that vectors in `A` are
uniformly distributed.
"""
function spacing(A::Array{Array{Float64,1},1})
    n = length(A)

    d = zeros(n)
    for i = 1:n
    	@inbounds d[i] = minimum( [ sum( abs.(A[i] - A[j]) ) for j in deleteat!(collect(1:n), i) ] )
    end

    d_hat = mean(d)
    s = (1.0 / ( n - 1.0 )) * sum( (d_hat .- d) .^ 2 )
    return sqrt(s)
    
end


spacing(A::Array{xFgh_indiv}) = spacing( fval.(A) )
spacing(A::State{xFgh_indiv}) = spacing(A.population)
spacing(A::Matrix) = spacing([A[i,:] for i in 1:size(A,1)])

"""
	covering(A, B)
Computes the covering indicator (percentage of vectors in B that are dominated by vectors in A)
from two sets with non-dominated solutions.

A and B with size (n, m) where n is number of samples and m is the vector dimension.

Note that `covering(A, B) == 1` means that all solutions in B are dominated by those
in A. Moreover, `covering(A, B) != covering(B, A)` in general.

If `A::State` and `B::State`, the computes `covering(A.population, B.population)` after
ignoring dominated solutions in each set.
"""
function covering(A::Array{Vector{T}}, B::Array{Vector{T}}) where T <: Real
    n_a = length(A)
    n_b = length(B)
    s = 0.0
    for i = 1:n_b
        for j = 1:n_a
			# `compare` returns 1 if argument 1 dominates argument 2
            if compare(A[j], B[i]) == 1
                s += 1.0
                break
            end    
        end
        
    end

    return s / n_b
    
end

covering(A::Vector{xFgh_indiv}, B::Vector{xFgh_indiv}) = covering(fval.(A), fval.(B))

function covering(A::State{xFgh_indiv}, B::State{xFgh_indiv})
	A_non_dominated = get_non_dominated_solutions(A.population)
	B_non_dominated = get_non_dominated_solutions(B.population)
	
	length(A_non_dominated) != length(A.population) &&
	@warn "Some solutions in B dominate each other (solutions ignored)."


	length(B_non_dominated) != length(B.population) &&
	@warn "Some solutions in B dominate each other (solutions ignored)."

	covering(A_non_dominated, B_non_dominated)
end


function covering(A::Matrix, B::Matrix)
	A_arr = [ A[i,:] for i in 1:size(A,1) ]
	B_arr = [ B[i,:] for i in 1:size(B,1) ]
	covering(A_arr, B_arr)
end


end
