module Metaheuristics

export WOA, CGSA, SA
export DE, ABC, PSO

# v2 items
export optimize, ECA, Options, State, Information, Problem
export MOEAD_DE

# v2.1 items
export convergence, minimizer, minimum, positions, fvals, nfes, get_position, fval, NSGA2
export hval, gval, hvals, gvals

export NSGA3, gen_ref_dirs

export SMS_EMOA, SPEA2, MCCGA, CCMO
export BRKGA

export PerformanceIndicators, pareto_front, nadir, ideal
export termination_status_message

export GA, RandomInBounds, RandomBinary, RandomPermutation
export TournamentSelection, RouletteWheelSelection
export UniformCrossover, OrderCrossover, SBX
export BitFlipMutation, SlightMutation, PolynomialMutation
export GenerationalReplacement, ElitistReplacement, RankAndCrowding
export sample, LatinHypercubeSampling, Grid 
export εDE, Restart
export optimize!
export set_user_solutions!
export boxconstraints
export SHADE

include("externals.jl")


include("core/abstracts.jl")
include("core/information.jl")
include("core/options.jl")
include("core/state.jl")
include("core/structures.jl")


include("solutions/individual.jl")
include("solutions/constrained.jl")
include("solutions/display.jl")

include("operators/operators.jl")

include("termination/termination.jl")

include("common/utils.jl")
include("common/multi-objective-functions.jl")
include("common/repair.jl")

include("common/compare.jl")
include("common/non-dominated-sorting.jl")
include("common/set_user_solutions.jl")



include("TestProblems/TestProblems.jl")
include("PerformanceIndicators/PerformanceIndicators.jl")


include("algorithms/algorithm.jl")
include("parameters/parameters.jl")

include("optimize/optimize.jl")

#######################################################
#                                                     #
#           A L G O R I T H M S                       #
#                                                     #
#######################################################

# template file
include("algorithms/template.jl")

include("algorithms/singleobjective/singleobjective.jl")
include("algorithms/multiobjective/multiobjective.jl")
include("algorithms/combinatorial/combinatorial.jl")

include("algorithms/stop_criteria.jl")
include("DecisionMaking/DecisionMaking.jl")
include("precompile/precompile.jl")

#######################################################
#                                                     #
#           D E P R E C A T E D                       #
#                                                     #
#######################################################
# nothing to deprecate yet

end # module
