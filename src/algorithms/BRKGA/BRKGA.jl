#struct BRKGA <: AbstractGA end

function BRKGA(;
        num_elites = 20,
        num_mutants = 10,
        num_offsprings = 70,
        N = num_elites + num_mutants + num_offsprings,
        bias = 0.7,
        kargs...
    )

    initializer = RandomInBounds(;N)
    selection = BiasedSelection(;num_elites = num_offsprings)
    crossover = BinomialCrossover(p = bias, n_offsprings = 1)
    mutation  = InsertRandomMutation(num_mutants)
    environmental_selection = ElitistReplacement()

    GA(;initializer,
       selection,
       crossover,
       mutation,
       environmental_selection,
       kargs...
      )
end

#=
function get_parameters(
        f,
        search_space::BoxConstrainedSpace,
        ::Type{T}
    ) where T <: BRKGA
    BRKGA()
end
=#
