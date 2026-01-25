# All functions are available in Metaheuristics namespace
import Random: AbstractRNG, rand, randn

"""
    LSRTDE(;N=0, information=Information(), options=Options())

Parameters for L-SRTDE (Linear population size reduction Success Rate-based adaptive Differential Evolution).

# Example

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], LSRTDE())
+=========== RESULT ==========+
  iteration: 1000
    minimum: 0
  minimizer: [0.0, 0.0, 0.0]
    f calls: 30000
 total time: 0.0437 s
+============================+
```
"""
mutable struct LSRTDE <: AbstractDifferentialEvolution
    N::Int
    N_front::Int
    N_front_max::Int
    MemorySize::Int
    MemoryCr::Vector{Float64}
    MemoryIter::Int
    SuccessRate::Float64
    F_sigma::Float64
    CR_sigma::Float64
    SuccessFilled::Int
    tempSuccessCr::Vector{Float64}
    FitDelta::Vector{Float64}
    front_population::Vector  # Front population (subset)
    front_fitness::Vector{Float64}  # Fitness of front population
    PFIndex::Int  # Index for cycling through front population
end

function LSRTDE(;
    N::Int = 0,
    kargs...
)
    MemorySize = 5
    parameters = LSRTDE(
        N,
        0,  # N_front (will be set in initialize!)
        0,  # N_front_max (will be set in initialize!)
        MemorySize,
        zeros(Float64, MemorySize),  # MemoryCr (will be initialized in initialize!)
        0,  # MemoryIter (0-indexed to match C++)
        0.5,  # SuccessRate
        0.02,  # F_sigma
        0.05,  # CR_sigma
        0,  # SuccessFilled
        Float64[],  # tempSuccessCr (will be resized in initialize!)
        Float64[],  # FitDelta (will be resized in initialize!)
        [],  # front_population
        Float64[],  # front_fitness
        0  # PFIndex (0-indexed to match C++, but will use 1-based indexing in Julia)
    )

    Algorithm(parameters; kargs...)
end

function initialize!(
    status,
    parameters::LSRTDE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    D = getdim(problem)

    if parameters.N <= 5
        parameters.N = 20 * D
    end

    # Initialize front population size
    parameters.N_front = parameters.N
    parameters.N_front_max = parameters.N

    # Initialize CR memory with ones
    parameters.MemoryCr = ones(Float64, parameters.MemorySize)
    parameters.MemoryIter = 0  # 0-indexed to match C++
    parameters.SuccessRate = 0.5
    parameters.SuccessFilled = 0

    # Resize temporary arrays
    parameters.tempSuccessCr = zeros(Float64, parameters.N * 2)
    parameters.FitDelta = zeros(Float64, parameters.N * 2)

    parameters.PFIndex = 0  # 0-indexed conceptually, but will use 1-based for Julia arrays

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000 * D
        options.debug &&
        @warn("f_calls_limit increased to $(options.f_calls_limit)")
    end

    if options.iterations == 0
        options.iterations = div(options.f_calls_limit, parameters.N) + 1
    end

    # Generate initial state
    state = gen_initial_state(problem, parameters, information, options, status)

    # Initialize front population as sorted copy of initial population
    front_pop = copy(state.population)
    sort!(front_pop, lt=(a, b) -> is_better(a, b, parameters))
    parameters.front_population = front_pop
    parameters.front_fitness = [fval(sol) for sol in front_pop]

    return state
end

"""
    mean_weighted_lehmers(Vector, Weights)

Calculate weighted Lehmer mean: SumSquare / Sum where:
- SumSquare = sum(Weights[i] * Vector[i]^2)
- Sum = sum(Weights[i] * Vector[i])
Returns 1.0 if Sum is too small.
"""
function mean_weighted_lehmers(Vector::Vector{Float64}, Weights::Vector{Float64}, n::Int)
    if n == 0
        return 1.0
    end

    # Normalize weights
    SumWeight = sum(Weights[1:n])
    if SumWeight == 0.0
        return 1.0
    end

    normalized_weights = Weights[1:n] ./ SumWeight

    # Calculate weighted sums
    SumSquare = sum(normalized_weights[i] * Vector[i] * Vector[i] for i in 1:n)
    Sum = sum(normalized_weights[i] * Vector[i] for i in 1:n)

    if abs(Sum) > 1e-8
        return SumSquare / Sum
    else
        return 1.0
    end
end

"""
    remove_worst!(front_population, front_fitness, new_size)

Remove worst individuals from front until size matches new_size.
"""
function remove_worst!(
    parameters::LSRTDE,
    new_size::Int
)
    current_size = length(parameters.front_population)
    if current_size <= new_size
        return
    end

    points_to_remove = current_size - new_size

    for _ in 1:points_to_remove
        # Find worst individual
        worst_idx = argmax(parameters.front_fitness)
        
        # Remove it
        deleteat!(parameters.front_population, worst_idx)
        deleteat!(parameters.front_fitness, worst_idx)
    end
end

function update_memory_cr!(parameters::LSRTDE)
    if parameters.SuccessFilled > 0
        new_cr = mean_weighted_lehmers(
            parameters.tempSuccessCr,
            parameters.FitDelta,
            parameters.SuccessFilled
        )
        # Use 1-based indexing for Julia arrays (MemoryIter is 0-indexed conceptually)
        parameters.MemoryCr[parameters.MemoryIter + 1] = 0.5 * (new_cr + parameters.MemoryCr[parameters.MemoryIter + 1])
        parameters.MemoryIter = (parameters.MemoryIter + 1) % parameters.MemorySize
    end
end

function update_state!(
    status,
    parameters::LSRTDE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    population = status.population
    N = parameters.N
    N_front = parameters.N_front
    D = getdim(problem)
    rng = options.rng

    # Track current population size (can grow with successful mutations)
    NIndsCurrent = length(population)
    
    # Reset success tracking
    parameters.SuccessFilled = 0

    # Calculate adaptive F mean
    meanF = 0.4 + tanh(parameters.SuccessRate * 5) * 0.25

    # Sort main population indices (only first N_front are considered for selection)
    main_fitness = [fval(sol) for sol in population]
    main_sorted_indices = sortperm(main_fitness[1:min(N_front, length(population))], rev=false)
    if length(main_sorted_indices) < N_front
        # Pad with indices if needed
        append!(main_sorted_indices, collect((length(main_sorted_indices)+1):N_front))
    end
    
    # Sort front population and update it to maintain sorted order
    front_sorted_indices = sortperm(parameters.front_fitness, rev=false)
    parameters.front_population = [parameters.front_population[i] for i in front_sorted_indices]
    parameters.front_fitness = [parameters.front_fitness[i] for i in front_sorted_indices]

    # Create discrete distribution weights for front selection
    front_weights = [exp(-i / N_front * 3) for i in 0:(N_front-1)]
    front_weights_sum = sum(front_weights)
    front_weights_normalized = front_weights ./ front_weights_sum

    # Calculate psizeval (top percentage of main population)
    psizeval = max(2, Int(floor(N_front * 0.7 * exp(-parameters.SuccessRate * 7))))

    # Generate trial vectors for each individual in front
    for i in 1:N_front
        # Select TheChosenOne randomly from front
        TheChosenOne = rand(rng, 1:N_front)
        x_chosen = get_position(parameters.front_population[TheChosenOne])

        # Select prand from top psizeval of sorted main population
        # C++ ensures prand != TheChosenOne
        prand_pos = rand(rng, 1:min(psizeval, length(main_sorted_indices)))
        prand_idx = main_sorted_indices[prand_pos]
        while prand_idx > length(population) || prand_idx == TheChosenOne
            prand_pos = rand(rng, 1:min(psizeval, length(main_sorted_indices)))
            prand_idx = main_sorted_indices[prand_pos]
        end
        x_prand = get_position(population[prand_idx])

        # Select Rand1 from front using discrete distribution
        # C++: Rand1 = Indices2[ComponentSelectorFront(...)] where Indices2[i] = i (sorted front indices)
        # C++ ensures Rand1 != prand where both are in range 0..NIndsFront-1
        # In C++, prand comes from Indices[IntRandom(psizeval)] which is an index into Popul
        # and Rand1 comes from Indices2[...] which is an index into PopulFront
        # The constraint ensures they don't select from the same position in their sorted arrays
        r = rand(rng)
        cumsum = 0.0
        Rand1_idx = 1
        for j in 1:N_front
            cumsum += front_weights_normalized[j]
            if r <= cumsum
                Rand1_idx = j
                break
            end
        end
        # Ensure Rand1 is different from prand position (matching C++ constraint)
        # prand_pos is 1..psizeval, Rand1_idx is 1..N_front
        # We ensure they're different positions when both are valid
        while Rand1_idx == prand_pos && N_front > 1 && prand_pos <= N_front
            r = rand(rng)
            cumsum = 0.0
            for j in 1:N_front
                cumsum += front_weights_normalized[j]
                if r <= cumsum
                    Rand1_idx = j
                    break
                end
            end
        end
        x_rand1 = get_position(parameters.front_population[Rand1_idx])

        # Select Rand2 from main population (using sorted indices)
        # C++ ensures Rand2 != prand && Rand2 != Rand1
        # Rand2 comes from Indices[IntRandom(NIndsFront)] (same array as prand)
        # So we compare indices directly
        Rand2_pos = rand(rng, 1:min(N_front, length(main_sorted_indices)))
        Rand2_idx = main_sorted_indices[Rand2_pos]
        while (Rand2_idx > length(population)) || Rand2_idx == prand_idx || Rand2_pos == Rand1_idx
            Rand2_pos = rand(rng, 1:min(N_front, length(main_sorted_indices)))
            Rand2_idx = main_sorted_indices[Rand2_pos]
        end
        x_rand2 = get_position(population[Rand2_idx])

        # Generate F from normal distribution
        F = randn(rng) * parameters.F_sigma + meanF
        while F < 0.0 || F > 1.0
            F = randn(rng) * parameters.F_sigma + meanF
        end
        F = clamp(F, 0.0, 1.0)

        # Generate CR from normal distribution using memory
        # C++ uses 0-based indexing: IntRandom(MemorySize) gives 0..MemorySize-1
        MemoryCurrentIndex = rand(rng, 0:parameters.MemorySize-1)
        Cr = randn(rng) * parameters.CR_sigma + parameters.MemoryCr[MemoryCurrentIndex + 1]  # +1 for 1-based Julia arrays
        Cr = clamp(Cr, 0.0, 1.0)

        # Apply L-SRTDE mutation strategy
        Trial = copy(x_chosen)
        ActualCr = 0.0
        j_rand = rand(rng, 1:D)

        for j in 1:D
            if rand(rng) < Cr || j == j_rand
                Trial[j] = x_chosen[j] + F * (x_prand[j] - x_chosen[j]) + F * (x_rand1[j] - x_rand2[j])
                ActualCr += 1.0
            end
        end
        ActualCr = ActualCr / D

        # Repair boundary violations (use reset_to_violated_bounds for better convergence)
        replace_with_random_in_bounds!(Trial, problem.search_space)

        # Evaluate trial solution
        trial_sol = create_solution(Trial, problem, ε=options.h_tol)
        trial_fit = fval(trial_sol)
        chosen_fit = parameters.front_fitness[TheChosenOne]

        # If better, replace and store for memory update
        if is_better(trial_sol, parameters.front_population[TheChosenOne], parameters)
            # Add to main population (growing it)
            push!(status.population, trial_sol)
            NIndsCurrent += 1

            # Replace using PFIndex cycling (matching C++ behavior)
            # PFIndex is 0-indexed conceptually, but we use 1-based for Julia arrays
            pf_idx = parameters.PFIndex + 1
            parameters.front_population[pf_idx] = trial_sol
            parameters.front_fitness[pf_idx] = trial_fit

            # Store successful CR and fitness delta
            if parameters.SuccessFilled < length(parameters.tempSuccessCr)
                parameters.tempSuccessCr[parameters.SuccessFilled + 1] = ActualCr
                parameters.FitDelta[parameters.SuccessFilled + 1] = abs(chosen_fit - trial_fit)
                parameters.SuccessFilled += 1
            end

            # Update PFIndex (cycling through front population)
            parameters.PFIndex = (parameters.PFIndex + 1) % N_front

            # Update best solution
            if is_better(trial_sol, status.best_sol, parameters)
                status.best_sol = trial_sol
            end
        end
    end

    # Update success rate
    parameters.SuccessRate = parameters.SuccessFilled / N_front

    # Calculate new front size (linear reduction)
    # C++ formula: int(double(4-NIndsFrontMax)/double(MaxFEval)*NFEval + NIndsFrontMax)
    # Note: C++ doesn't clamp to 4, but we'll keep the clamp for safety
    MaxFEval = options.f_calls_limit
    NFEval = status.f_calls
    newN_front = Int(floor((4.0 - parameters.N_front_max) / MaxFEval * NFEval + parameters.N_front_max))
    newN_front = max(4, newN_front)  # Clamp to at least 4 (C++ doesn't clamp, but this is safer)

    # Remove worst individuals from front
    remove_worst!(parameters, newN_front)
    parameters.N_front = length(parameters.front_population)

    # Update CR memory
    update_memory_cr!(parameters)

    # Sort and trim population to maintain size
    # Sort all individuals in population
    sort!(status.population, lt=(a, b) -> is_better(a, b, parameters))
    
    # Keep only best N_front
    if length(status.population) > parameters.N_front
        status.population = status.population[1:parameters.N_front]
    end
    
    # Update front to match sorted population (front should be best N_front)
    parameters.front_population = copy(status.population)
    parameters.front_fitness = [fval(sol) for sol in status.population]
    parameters.N_front = length(parameters.front_population)

    # Reset PFIndex for next generation (C++ sets it to 0 at start of each generation cycle)
    parameters.PFIndex = 0

    # Reset SuccessFilled for next generation
    parameters.SuccessFilled = 0

    # Update best solution
    status.best_sol = get_best(status.population)
end

function final_stage!(
    status,
    parameters::LSRTDE,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
)
    status.final_time = time()
    status.best_sol = get_best(status.population)
end

is_better(a, b, parameters::LSRTDE) = is_better(a, b)
