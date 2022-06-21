
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
	igd_plus(front, true_pareto_front; p = 1)

Returns the Inverted Generational Distance Plus.

### Parameters

`front` and `true_pareto_front` can be:

- `N×m` matrix where `N` is the number of points and `m` is the number of objectives. 
- `State`
- `Array{xFgh_indiv}` (usually `State.population`)

"""
igd_plus(front, true_pareto_front; p = 1) = generational_distance(front, true_pareto_front, inverted=true, p=p, plus=true)

