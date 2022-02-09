function update_state!(
        status::State{xfgh_indiv},
        parameters::ECA,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    I = randperm(parameters.N)
    N = parameters.N
    population = status.population

    ε = 0.0

    if parameters.ε > 0.0
        p = status.f_calls / options.f_calls_limit
        ε = parameters.ε * (1-p)^4
    end

    feasible_solutions = findall( s->s.is_feasible, population )
    weights = compute_weights(population)

    X_next = zeros(parameters.N, size(problem.bounds, 2))

    # For each elements in Population
    for i = 1:parameters.N
        # current
        x = population[i].x

        # generate U masses
        U_ids = getU_ids(parameters.K, I, i, parameters.N, feasible_solutions)
        U = population[U_ids]

        # generate center of mass
        mass = weights[U_ids]
        c = center(U, mass)
        u_worst = argmin(mass)
        u_best = argmax(mass)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x
        v = U[u_best].x

        # current-to-center/bin
        y = x .+ η .* (c .- u)
        
        mask = rand(length(y)) .< 1.0 / length(y)
        y[mask] = v[mask]

        evo_boundary_repairer!(y, c, problem.bounds)
        X_next[i,:] = y
    end

    for sol in create_solutions(X_next, problem, ε=ε)
        if is_better(sol, status.best_sol)
            status.best_sol = sol
        end
        # stored but not used until replacement step
        push!(status.population, sol)
        status.stop = stop_check(status, information, options)
        status.stop && break
    end
    
    if length(status.population) == N
        # no resize population
        return
    end

    # replacement step
    sort!(status.population, lt = is_better, alg=QuickSort)
    deleteat!(status.population, N+1:length(status.population))
end

function compute_weights(population)
    fs = fvals(population)
    fx_max = maximum( abs.(fs) )
    vios = map(sum_violations, population)
    fs += 2.0*fx_max*vios
    γ = maximum( abs.(fs) )
    2.0γ .- fs
end

function final_stage!(
        status::State{xfgh_indiv},
        parameters::ECA,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    status.final_time = time()

end

