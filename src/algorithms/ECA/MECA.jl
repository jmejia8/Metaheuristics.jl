mutable struct MECA <: AbstractParameters
    η_max::Float64
    K::Int
    N::Int
    ε::Float64
    λ::Array{Vector{Float64}}
    z::Vector{Float64}
    nparitions::Int
    nobjectives
end

"""
    MECA(;
        η_max = 2.0,
        K = 7,
        N = 0,
        ε = 0.0,
        information = Information(),
        options = Options()
    )

Parameters for the metaheuristic ECA: step-size `η_max`,`K` is number of vectors to
generate the center of mass, `N` is the population size.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], MECA())

+=========== RESULT ==========+
| Iter.: 1021
| f(x) = 1.68681e-163
| solution.x = [2.5517634463667404e-82, -2.9182760041942484e-82, -1.3565584801935802e-82]
| f calls: 21454
| Total time: 0.0894 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], ECA(N = 10, η_max = 1.0, K = 3))
+=========== RESULT ==========+
| Iter.: 1506
| f(x) = 0.000172391
| solution.x = [-6.340714627875324e-5, -0.004127226953894587, 0.012464071313908906]
| f calls: 15069
| Total time: 0.0531 s
+============================+
```

"""
function MECA(;
    η_max::Float64 = 2.0,
    K::Int = 7,
    N::Int = 90,
    ε::Float64 = 0.5,
    z::Vector{Float64} = zeros(0),
    λ = Array{Vector{Float64}}[],
    nparitions = 12,
    information = Information(),
    options = Options(),
)

    parameters = MECA(
        η_max,
        K,
        N,
        ε, λ, z, nparitions, nparitions
    )
    alg = Algorithm(
        parameters,
        information = information,
        options = options,
    )


    alg

end


function update_state!(
    status::State,
    parameters::MECA,
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
    nrefs = length(parameters.λ)

    fs = map( i -> g(status.population[i].f, parameters.λ[ 1 + (i-1) % nrefs], parameters.z), 1:parameters.N )
    gs = copy(fs)

    fx_max = maximum( abs.(fs) )

    vios = map(sum_violations, population)
    fs += 2.0*fx_max*vios
    γ = maximum( abs.(fs) )
    weights = 2.0γ .- fs

    new_solutions = xFgh_indiv[]
    used = zeros(Int, length(parameters.λ))
    tmp_counter = 0

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

        sol = generateChild(y, problem.f(y), ε=options.h_tol)
        update_reference_point!(parameters.z, sol)
        status.f_calls += 1

        # stored but not used until replacement step
        j = i
        if j > nrefs
            j = rand(1:nrefs)
        end
        
        gsv = g(sol.f, parameters.λ[j], parameters.z)
        gs_old = g(status.population[i].f, parameters.λ[j], parameters.z)

        if gsv < gs_old
            used[j] += 1
            population[i] = sol
            tmp_counter += 1
        else
            push!(new_solutions, sol)
        end

        status.stop = stop_check(status, information, options)
        status.stop && break


        if tmp_counter >= parameters.N
            @info "Cool..."
            break
        end
        
    end

    if isempty(new_solutions)
        @info "wow"
        return
    end
    
    
    # save best solutions according to g and λ
    tmp_counter = 0 # count saved solutions per iteration
    gs_new = copy(gs)
    J = sortperm(used)[1:length(parameters.λ)]
    for j = J
        if used[j] > 0
            break
        end
        for i = 1:length(new_solutions)
            
            gsv = g(new_solutions[i].f, parameters.λ[j], parameters.z)
            gs_old = g(status.population[j].f, parameters.λ[j], parameters.z)
            if gsv < gs_old
                gs[j] = gsv
                population[j] = new_solutions[i]
                used[j] += 1
                tmp_counter+=1
                break
            end
        end
    end

    @show sum(used)/parameters.N
    @show used
    

    
end


function initialize!(
    status::State,
    parameters::MECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    D = size(problem.bounds, 2)


    if parameters.N <= 0
        error("Increase the population size (e.g. MECA(;N=90))")
    end

    if options.iterations == 0
        options.iterations = 500
    end

    if options.f_calls_limit == 0
        options.f_calls_limit = options.iterations * parameters.N + 1
    end


    initialize!(problem, nothing, parameters, status, information, options)
    parameters.nobjectives = length(status.population[1].f)
    parameters.z = fill(Inf, parameters.nobjectives)

    update_reference_point!(parameters.z, status.population)




    if isempty(parameters.λ)
        options.debug && @info "Generating reference points..."
        parameters.λ = gen_weights(parameters.nobjectives, parameters.nparitions)
    end

end

function final_stage!(
    status::State,
    parameters::MECA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    status.final_time = time()

    # compute Pareto front if it is a multiobjective problem
    if typeof(status.population[1].f) <: Array
        options.debug && @info "Computing Pareto front..."
        status.best_sol = get_pareto_front(status.population, is_better_eca)
    end
end

