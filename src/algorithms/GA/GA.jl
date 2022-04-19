abstract type AbstractGA <: AbstractParameters end

mutable struct GA{T1, T2, T3, T4, T5} <: AbstractGA
    initializer::T1
    selection::T2
    crossover::T3
    mutation::T4
    environmental_selection::T5
end

function GA(;
        N = 100,
        p_mutation  = 1e-5,
        p_crossover = 0.5,
        initializer = RandomInBounds(;N),
        selection   = TournamentSelection(;N),
        crossover   = UniformCrossover(p = p_crossover),
        mutation    = BitFlipMutation(p = p_mutation),
        environmental_selection = ElitistReplacement(),
        options     = Options(),
        information = Information()
    )
    
    parameters = GA(initializer,
                    selection,
                    crossover,
                    mutation,
                    environmental_selection
                   )

    Algorithm(
        parameters,
        information = information,
        options = options
    )
end


function initialize!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    @assert parameters.initializer.N > 0

    if options.iterations == 0
        options.iterations = 500
        options.f_calls_limit = 2options.iterations * parameters.initializer.N
    end

    initialize!(status,parameters.initializer,problem,information,options,args...;kargs...)
end


function update_state!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    population = status.population

    # parent selection
    parent_mask = selection(population, parameters.selection)
    # variation operators
    Q = crossover(population[parent_mask], parameters.crossover)
    mutation!(Q, parameters.mutation)
    # evaluate offspring
    offsprings = create_solutions(Q, problem)
    cb = get_best(status.population)
    # select solutions (survivals are saved in population)
    environmental_selection!(population, offsprings, parameters.environmental_selection)

    if is_better(cb, status.best_sol)
        # save elite solution
        status.best_sol = cb
    end
    
end



function final_stage!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    status.final_time = time()
    return
end

function stop_criteria!(
        status,
        parameters::AbstractGA,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    return
end

