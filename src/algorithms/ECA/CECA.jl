function center_ceca(U, max_fitness, ε)
    mass = getMass(U, max_fitness, ε)
    return center(U, mass),
    argworst(U), # worst
    argbest(U)  # best
end



function update_state!(
    status::State{xfgh_indiv},
    parameters::ECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    K = parameters.K
    I = randperm(parameters.N)
    N = parameters.N
    population = status.population
    D = size(problem.bounds, 2)


    a = problem.bounds[1, :]
    b = problem.bounds[2, :]

    ε = 0.0

    p = status.f_calls / options.f_calls_limit

    if parameters.ε > 0
        ε = parameters.ε * (p - 1)^4
    end


    feasible_solutions = findall( s->s.is_feasible, status.population )

    fs = fvals(status)
    fx_max = maximum( abs.(fs) )

    vios = map(sum_violations, population)
    fs += 2.0*fx_max*vios
    γ = maximum( abs.(fs) )
    weights = 2.0γ .- fs

    # For each elements in Population
    for i = 1:parameters.N

        # current
        x = status.population[i].x

        # generate U masses
        #U = getU(status.population, parameters.K, I, i, parameters.N, feasible_solutions)

        U_ids = getU_ids(parameters.K, I, i, parameters.N, feasible_solutions)
        U = population[U_ids]

        # generate center of mass
        mass = weights[U_ids]
        c = center(U, mass)
        u_worst = argmin(mass)

        # stepsize
        η = parameters.η_max * rand()

        # u: worst element in U
        u = U[u_worst].x

        # current-to-center/bin
        y = x .+ η .* (c .- u)

        evo_boundary_repairer!(y, c, problem.bounds)

        sol = create_solution(y, problem, ε=options.h_tol)
        status.f_calls += 1


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
    sort!(status.population, lt = is_better)
    deleteat!(status.population, N+1:length(status.population))
    


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

