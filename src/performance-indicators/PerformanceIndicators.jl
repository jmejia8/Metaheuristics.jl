module PerformanceIndicators

import ..State, ..fvals, ..norm

"""
	generational_distance(front, true_pareto_front; p = 1, inverted = false, plus=false)

Returns the Generational Distance.

### Parameters

- `front` is `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
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


"""
	gd(front, true_pareto_front; p = 1)

Returns the Generational Distance.

### Parameters

- `front` is `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `true_pareto_front` is a `M×m` matrix.
"""
gd(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, p=p)
gd(front::State, true_pareto_front::State; p = 1) = gd(fvals(front), fvals(true_pareto_front), p=p)


"""
	igd(front, true_pareto_front; p = 1)

Returns the Inverted Generational Distance.

### Parameters

- `front` is `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `true_pareto_front` is a `M×m` matrix.
"""
igd(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, inverted=true, p=p)
igd(front::State, true_pareto_front::State; p = 1) = igd(fvals(front), fvals(true_pareto_front), p=p)

"""
	gd_plus(front, true_pareto_front; p = 1)

Returns the Generational Distance Plus.

### Parameters

- `front` is `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `true_pareto_front` is a `M×m` matrix.
"""
gd_plus(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, p=p, plus=true)
gd_plus(front::State, true_pareto_front::State; p = 1) = igd_plus(fvals(front), fvals(true_pareto_front), p=p)

"""
	igd_plus(front, true_pareto_front; p = 1)

Returns the Inverted Generational Distance Plus.

### Parameters

- `front` is `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `true_pareto_front` is a `M×m` matrix.
"""
igd_plus(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, inverted=true, p=p, plus=true)
igd_plus(front::State, true_pareto_front::State; p = 1) = igd_plus(fvals(front), fvals(true_pareto_front), p=p)


end
