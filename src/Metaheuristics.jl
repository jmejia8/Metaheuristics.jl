__precompile__()
module Metaheuristics

export WOA, GOA, GSA, CGSA, SA, CMAES_AEP
export CFO, DE, ABC, PSO

# v2 items
export optimize, ECA, Options, State, Information, Engine, Problem
export MOEAD_DE

# v2.1 items
export convergence, minimizer, minimum, positions, fvals, nfes, get_position, fval, NSGA2
export hval, gval, hvals, gvals

export NSGA3, gen_ref_dirs

export SMS_EMOA, SPEA2, MCCGA

export PerformanceIndicators, pareto_front, nadir, ideal
export termination_status_message

export GA, RandomInBounds, RandomBinary, RandomPermutation
export TournamentSelection, RouletteWheelSelection
export UniformCrossover, OrderCrossover, SBX
export BitFlipMutation, SlightMutation, PolynomialMutation
export GenerationalReplacement, ElitistReplacement

include("externals.jl")


include("core/abstracts.jl")
include("core/stop_status_codes.jl")
include("core/information.jl")
include("core/options.jl")
include("core/state.jl")
include("core/structures.jl")


include("solutions/individual.jl")
include("solutions/constrained.jl")
include("solutions/display.jl")

include("operators/operators.jl")

include("common/utils.jl")
include("common/multi-objective-functions.jl")
include("common/repair.jl")

include("common/stop.jl")
include("common/compare.jl")
include("common/non-dominated-sorting.jl")


include("TestProblems/TestProblems.jl")
include("PerformanceIndicators/PerformanceIndicators.jl")


include("algorithms/algorithm.jl")

include("optimize.jl")

#######################################################
#                                                     #
#           A L G O R I T H M S                       #
#                                                     #
#######################################################

# template file
include("algorithms/template.jl")

include("algorithms/GA/GA.jl")

# ECA algorithm
include("algorithms/ECA/ECA.jl")
include("algorithms/ECA/CECA.jl")

# Differential algorithm
include("algorithms/DE/DE.jl")

# PSO algorithm
include("algorithms/PSO/PSO.jl")

# The whale optimization algorithm
# S Mirjalili, A Lewis - Advances in Engineering Software, 2016
include("algorithms/WOA/WOA.jl")

# Mirjalili, Seyedali, and Amir H. Gandomi.
# "Chaotic gravitational constants for the gravitational search algorithm."
# Applied Soft Computing 53 (2017): 407-419.
include("algorithms/CGSA/CGSA.jl")

# SA: Simulated Annealing
# Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
# Simulated Annealing. _Science, 220_, 671-680.
include("algorithms/SA/SA.jl")


# Aritifical Bee colony
include("algorithms/ABC/ABC.jl")

# MOEA decomposition based using Differential Evolution
include("algorithms/MOEAD_DE/MOEAD_DE.jl")

# Non-dominate sorting Genetic Algorithm
include("algorithms/NSGA2/NSGA2.jl")
include("algorithms/NSGA3/NSGA3.jl")
include("algorithms/SMS_EMOA/SMS_EMOA.jl")
include("algorithms/SPEA2/SPEA2.jl")

# genetic algorithm
include("algorithms/MCCGA/MCCGA.jl")


include("algorithms/stop_criteria.jl")

include("DecisionMaking/DecisionMaking.jl")

#######################################################
#                                                     #
#           D E P R E C A T E D                       #
#                                                     #
#######################################################
# nothing to deprecate yet

end # module
