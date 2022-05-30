struct ROIArchiving <: AbstractDecisionMakingMethod
    δ_w::Vector{Float64}
    δ_w_transformed::Vector{Float64}

    ROIArchiving(δ_w) = begin
        δ_w_transformed = δ_w / 2.0
        δ_w_transformed = 1.0 .- cos.(π/2*δ_w_transformed)
        new(δ_w,  δ_w_transformed)
    end    
end

function roiarchiving(population, w, parameters::ROIArchiving; verbose=true)

    if isempty(population) || isempty(w)
        # nothing to do
        return Int[]
    end

    feasible_sols = Metaheuristics.is_feasible.(population)
    if !any(feasible_sols)
        # we cannot work on infeasible solutions
        return Int[]
    end
    
    δ_w = parameters.δ_w_transformed

    size(w,2) != length(δ_w) && error("|weight_points| is different to |δ_w|")

    non_dominated = Metaheuristics.get_non_dominated_solutions_perm(population)
    isempty(non_dominated) && return Int[]

    population = population[non_dominated]

    fmin = ideal(population)
    fmax = nadir(population)

    g_roi = compute_rio_vio(population, fmin, fmax, w, δ_w)

    # return preferred solutions in Region of Interest
    idx = findall(<=(zero(eltype(g_roi))), g_roi)
    if !isempty(idx)
        return idx
    end

    verbose && @warn "Returning the feasible solution closest to the region of interest."
    return [argmin(g_roi)]

end

cosine_distance(x, y) = one(eltype(x)) - dot(x, y) / (norm(x) * norm(y))

function compute_rio_vio(population, fmin, fmax, w, δ_w)
    g = zeros(length(population))
    n = 1:size(w,1)
    Δ = (fmax - fmin)
    for (l,s) in enumerate(population)
        # filtering based on weight_points
        gx = minimum( i -> cosine_distance(view(w, i, :), (fval(s) - fmin) ./ Δ) - δ_w[i], n)
        g[l] = max(gx,0)
    end
    return g
end
