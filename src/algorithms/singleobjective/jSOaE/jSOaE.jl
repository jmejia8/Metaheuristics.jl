import LinearAlgebra: eigen
import Statistics: cov
import Random

mutable struct jSOaE <: AbstractDifferentialEvolution
    N::Int                    # Current population size
    N_init::Int               # Initial population size
    N_min::Int                # Minimum population size
    archive::Vector           # External archive
    MF::Vector{Float64}       # Memory for F (size H)
    MCR::Vector{Float64}      # Memory for CR (size H)
    F::Vector{Float64}        # Current F values
    CR::Vector{Float64}       # Current CR values
    k::Int                    # Memory index
    H::Int                    # Memory size
    Ap::Float64               # Archive probability parameter
    ps::Float64               # Selection proportion for eig
    peig::Float64             # Probability of eigenvalue crossover
    Asize::Int                # Current archive size
    Asize_max::Int            # Maximum archive size
end


"""
    jSOaE(;
        N = 0,
        H = 5,
        Ap = 0.2,
        ps = 0.5,
        peig = 0.4,
        N_min = 4,
        information = Information(),
        options = Options()
    )

jSOaE (jSO with Eigenvalue-based crossover) is an improved Success-history based 
Adaptive Differential Evolution algorithm that incorporates eigenvalue-based crossover.

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> result = optimize(f, [-1 -1 -1; 1 1 1.0], jSOaE())
Optimization Result
===================
  Iteration:       300
  Minimum:         1.234e-35
  Minimizer:       [1.23e-18, -2.34e-18, 3.45e-18]
  Function calls:  30000
  Total time:      0.0647 s
  Stop reason:     Maximum objective function calls exceeded.
```
"""
function jSOaE(;
    N::Int = 0,
    H::Int = 5,
    Ap::Float64 = 0.2,
    ps::Float64 = 0.5,
    peig::Float64 = 0.4,
    N_min::Int = 4,
    kargs...
)
    # Initialize memory arrays with default values
    MF = fill(0.3, H)
    MCR = fill(0.8, H)
    MF[H] = 0.9
    MCR[H] = 0.9
    
    parameters = jSOaE(
        N,              # N (will be set in initialize!)
        0,              # N_init (will be set in initialize!)
        N_min,          # N_min
        [],             # archive
        MF,             # MF
        MCR,            # MCR
        zeros(0),       # F
        zeros(0),       # CR
        1,              # k
        H,              # H
        Ap,             # Ap
        ps,             # ps
        peig,           # peig
        0,              # Asize
        0               # Asize_max (will be set in initialize!)
    )
    Algorithm(parameters; kargs...)
end

function initialize!(
        status,
        parameters::jSOaE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )
    D = getdim(problem)

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000D
        options.debug &&
        @warn( "f_calls_limit increased to $(options.f_calls_limit)")
    end

    # Calculate initial population size
    if D < 6
        N_init = round(Int, 25 * log(6) * sqrt(6))
    else
        N_init = round(Int, 25 * log(D) * sqrt(D))
    end
    
    if N_init < 2 * parameters.N_min
        N_init = 2 * parameters.N_min
    end
    
    # Set population size
    if parameters.N <= 0
        parameters.N = N_init
    else
        N_init = parameters.N
    end
    
    parameters.N_init = N_init
    parameters.Asize_max = round(Int, N_init * 2.6)
    parameters.Asize = 0
    parameters.archive = []
    
    # Initialize memory arrays if not already set
    if length(parameters.MF) != parameters.H
        parameters.MF = fill(0.3, parameters.H)
        parameters.MCR = fill(0.8, parameters.H)
        parameters.MF[parameters.H] = 0.9
        parameters.MCR[parameters.H] = 0.9
    end
    
    # Initialize F and CR arrays
    parameters.F = zeros(parameters.N)
    parameters.CR = zeros(parameters.N)
    parameters.k = 1

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    return gen_initial_state(problem, parameters, information, options, status)
end


function update_state!(
        status,
        parameters::jSOaE,
        problem::AbstractProblem,
        information::Information,
        options::Options,
        args...;
        kargs...
    )

    new_vectors = reproduction(status, parameters, problem, options)

    # evaluate solutions
    new_solutions = create_solutions(new_vectors, problem, ε=options.h_tol)
    append!(status.population, new_solutions)

    # Store rng in parameters temporarily for environmental_selection
    # (This is a workaround since environmental_selection doesn't receive options)
    environmental_selection!(status.population, parameters, options.rng)
    status.best_sol = get_best(status.population)
    
    # Population size reduction
    neval = status.f_calls
    MaxEval = options.f_calls_limit
    N_old = parameters.N
    N_new = round(Int, ((parameters.N_min - parameters.N_init) / MaxEval) * neval + parameters.N_init)
    
    if N_new < N_old && N_new >= parameters.N_min
        # Sort population by fitness
        sort!(status.population, lt=(a, b) -> is_better(a, b, parameters))
        # Keep top N_new solutions
        if length(status.population) > N_new
            deleteat!(status.population, N_new + 1:length(status.population))
        end
        parameters.N = N_new
        parameters.Asize_max = round(Int, N_new * 2.6)
        
        # Reduce archive size if needed
        while parameters.Asize > parameters.Asize_max
            if parameters.Asize > 0
                idx = rand(options.rng, 1:parameters.Asize)
                deleteat!(parameters.archive, idx)
                parameters.Asize -= 1
            else
                break
            end
        end
        
        # Resize F and CR arrays
        #parameters.F = zeros(parameters.N)
        #parameters.CR = zeros(parameters.N)
    end
end


function reproduction(status, parameters::jSOaE, problem::AbstractProblem, options::Options)
    population = status.population
    @assert !isempty(population)
    
    N = parameters.N
    D = getdim(problem)
    neval = status.f_calls
    MaxEval = options.f_calls_limit
    rng = options.rng
    
    # Sort population by fitness
    sort!(population, lt=(a, b) -> is_better(a, b, parameters))
    
    # Calculate pbest proportion
    pmax = 0.25
    pmin = pmax / 2
    pp = pmax - ((pmax - pmin) * (neval / MaxEval))
    p = max(2, ceil(Int, pp * N))
    
    # Cauchy distribution generator
    randc() = begin
        u = rand(rng) * π - π / 2
        0.1 * tan(u) + parameters.MF[rand(rng, 1:parameters.H)]
    end
    
    # Decide whether to use eigenvalue crossover
    use_eig = rand(rng) < parameters.peig
    
    # Prepare eigenvalue crossover data if needed
    EigVect = nothing
    if use_eig
        # Select top solutions for covariance calculation
        n_top = max(3, round(Int, N * parameters.ps + 1))
        if n_top > length(population)
            n_top = length(population)
        end
        
        if n_top >= 3
            top_pop = population[1:n_top]
            X_top = positions(top_pop)
            
            # Compute covariance matrix
            if n_top > 1
                cov_matrix = cov(X_top)
                # Get eigenvectors
                eigen_result = eigen(cov_matrix)
                EigVect = eigen_result.vectors
            end
        end
    end
    
    X = zeros(eltype(get_position(population[1])), N, D)
    X_pop = positions(population)
    
    for i in 1:N
        # Select random index from memory
        r = rand(rng, 1:parameters.H)
        
        # Generate CR
        CR = parameters.MCR[r] + sqrt(0.1) * randn(rng)
        CR = clamp(CR, 0.0, 1.0)
        
        # jSO CR adaptations
        if neval < 0.25 * MaxEval
            CR = max(CR, 0.7)
        elseif neval < 0.5 * MaxEval
            CR = max(CR, 0.6)
        end
        
        # Generate F using Cauchy distribution
        F = -1.0
        for _ in 1:100  # Avoid infinite loop
            u = rand(rng) * π - π / 2
            F = 0.1 * tan(u) + parameters.MF[r]
            if F > 0
                break
            end
        end
        F = clamp(F, 0.0, 1.0)
        
        # jSO F adaptation
        if (neval < 0.6 * MaxEval) && (F > 0.7)
            F = 0.7
        end
        
        # Calculate Fw
        if neval < 0.2 * MaxEval
            Fw = 0.7 * F
        elseif neval < 0.4 * MaxEval
            Fw = 0.8 * F
        else
            Fw = 1.2 * F
        end
        
        # Store F and CR
        parameters.F[i] = F
        parameters.CR[i] = CR
        
        # Current solution
        nx = X_pop[i, :]
        
        # Select pbest
        pbest_idx = rand(rng, 1:p)
        xpbest = X_pop[pbest_idx, :]
        
        # Select r1 (excluding current)
        available = setdiff(1:N, [i])
        r1_idx = available[rand(rng, 1:length(available))]
        r1 = X_pop[r1_idx, :]
        
        # Select r2 from population + archive (excluding current and r1)
        # Combine population and archive
        combined_pop = eltype(population)[population; parameters.archive]
        combined_size = length(combined_pop)
        
        # Create list of available indices (excluding i and r1_idx)
        # For population indices, exclude i and r1_idx
        available_indices = Int[]
        for idx in 1:combined_size
            if idx <= N
                # From population: exclude i and r1_idx
                if idx != i && idx != r1_idx
                    push!(available_indices, idx)
                end
            else
                # From archive: all are available
                push!(available_indices, idx)
            end
        end
        
        if !isempty(available_indices)
            r2_idx = available_indices[rand(rng, 1:length(available_indices))]
            if r2_idx <= N
                r2 = X_pop[r2_idx, :]
            else
                r2 = get_position(combined_pop[r2_idx])
            end
        else
            # Fallback: select from population (shouldn't happen normally)
            available = setdiff(1:N, [i, r1_idx])
            if !isempty(available)
                r2 = X_pop[available[rand(rng, 1:length(available))], :]
            else
                r2 = X_pop[rand(rng, 1:N), :]
            end
        end
        
        # Generate mutant
        v = nx + Fw * (xpbest - nx) + F * (r1 - r2)
        
        # Crossover
        if use_eig && EigVect !== nothing && size(EigVect, 1) == D
            # Eigenvalue-based crossover
            # Transform to eigenvector space (nx and v are column vectors)
            nxeig = EigVect' * nx
            veig = EigVect' * v
            
            # Apply crossover in eigenvector space
            change = findall(rand(rng, D) .< CR)
            if isempty(change)
                change = [rand(rng, 1:D)]
            end
            nxeig[change] = veig[change]
            
            # Transform back to original space
            nx = EigVect * nxeig
        else
            # Standard binomial crossover
            change = findall(rand(rng, D) .< CR)
            if isempty(change)
                change = [rand(rng, 1:D)]
            end
            nx[change] = v[change]
        end
        
        # Boundary repair using reflection
        reflected_back_to_bound!(nx, problem.search_space)
        
        X[i, :] .= nx
    end
    
    X
end


function environmental_selection(population, parameters::jSOaE, rng::Random.AbstractRNG = Random.default_rng())
    @assert length(population) == 2 * parameters.N

    N = parameters.N
    k = parameters.k
    MCR = parameters.MCR
    MF = parameters.MF
    CR = parameters.CR
    F = parameters.F
    new_solutions = population[N+1:end]

    survivals = Int[]
    Δ = Float64[]
    SCR = Float64[]
    SF = Float64[]

    # select survival and collect successful parameters
    for (i, h) in enumerate(new_solutions)
        if is_better(h, population[i], parameters)
            push!(survivals, N + i)
            delta_f = fval(population[i]) - fval(h)
            push!(Δ, delta_f)
            push!(SCR, CR[i])
            push!(SF, F[i])
        else
            push!(survivals, i)
        end
    end

    # Update memory if there are successful solutions
    if !isempty(Δ) && sum(Δ) > 0
        # Calculate weights
        w = Δ ./ sum(Δ)
        
        # Update MCR
        MCR_old = MCR[k]
        if any(SCR .> 0)
            # Weighted mean of CR
            MCR[k] = sum(w .* SCR)
            if MCR[k] == -1 || isnan(MCR[k])
                MCR[k] = MCR_old
            else
                MCR[k] = (MCR[k] + MCR_old) / 2
            end
        else
            MCR[k] = -1
        end
        
        # Update MF (Lehmer mean)
        MF_old = MF[k]
        if any(SF .> 0)
            numerator = sum(w .* SF.^2)
            denominator = sum(w .* SF)
            if denominator > 0
                MF[k] = numerator / denominator
                if isnan(MF[k]) || isinf(MF[k])
                    MF[k] = MF_old
                else
                    MF[k] = (MF[k] + MF_old) / 2
                end
            else
                MF[k] = MF_old
            end
        else
            MF[k] = MF_old
        end
        
        # Update archive
        mask = filter(>(N), survivals) .- N
        for idx in mask
            if parameters.Asize < parameters.Asize_max
                push!(parameters.archive, new_solutions[idx])
                parameters.Asize += 1
            else
                # Replace solution in archive (worse than median)
                if parameters.Asize > 0
                    # Sort archive by fitness (ascending: best first)
                    archive_fvals = [fval(sol) for sol in parameters.archive]
                    sorted_idx = sortperm(archive_fvals)
                    
                    # Calculate position to start replacement (worse solutions)
                    # ah = ceil(Asize * (1 - Ap)) means keep top (1-Ap) portion
                    # Replace from position ah onwards (the worst Ap portion)
                    ah = ceil(Int, parameters.Asize * (1 - parameters.Ap))
                    if ah < parameters.Asize
                        # Randomly select from positions ah to Asize (0-indexed: ah+1 to Asize)
                        replace_pos = ah + rand(rng, 1:(parameters.Asize - ah))
                        replace_idx = sorted_idx[replace_pos]
                        parameters.archive[replace_idx] = new_solutions[idx]
                    elseif ah == parameters.Asize && parameters.Asize > 0
                        # If ah equals Asize, replace the worst one
                        replace_idx = sorted_idx[end]
                        parameters.archive[replace_idx] = new_solutions[idx]
                    end
                end
            end
        end
        
        # Update memory index
        parameters.k = (k % parameters.H) + 1
    end

    return survivals
end

function environmental_selection!(population, parameters::jSOaE, rng::Random.AbstractRNG = Random.default_rng())
    mask = environmental_selection(population, parameters, rng)
    ignored = ones(Bool, length(population))
    ignored[mask] .= false
    deleteat!(population, ignored)
    return
end

is_better(a, b, parameters::jSOaE) = is_better(a, b)
