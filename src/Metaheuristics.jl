__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso, WOA, GOA, GSA, CGSA, SA, CMAES_AEP
export CFO, jso, DE, ABC, PSO

# v2 items
export optimize, ECA, Options, State, Information, Engine, Problem
export MOEAD_DE

# v2.1 items
export convergence, minimizer, minimum, positions, fvals, nfes, get_position, fval, NSGA2

include("externals.jl")


include("core/abstracts.jl")
include("core/information.jl")
include("core/options.jl")
include("core/state.jl")
include("core/structures.jl")


include("common/multi-objective-functions.jl")
include("common/repair.jl")
include("common/stop.jl")

include("solutions/individual.jl")
include("solutions/constrained.jl")
include("solutions/display.jl")
include("solutions/methods.jl")



include("algorithms/algorithm.jl")

include("optimize.jl")

#######################################################
#                                                     #
#           A L G O R I T H M S                       #
#                                                     #
#######################################################

# ECA algorithm
include("algorithms/ECA/ECA.jl")

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




#######################################################
#                                                     #
#           D E P R E C A T E D                       #
#                                                     #
#######################################################
include("deprecated/tools.jl")


# Grasshopper optimisation algorithm: Theory and application
# S Saremi, S Mirjalili, A Lewis - Advances in Engineering Software, 2017
include("deprecated/GOA.jl")

# GSA: a gravitational search algorithm
# E Rashedi, H Nezamabadi-Pour, S Saryazdi - Information sciences, 2009
include("deprecated/GSA.jl")

# Li, Zhenhua, and Qingfu Zhang.
# "An efficient rank-1 update for Cholesky CMA-ES using auxiliary evolution path."
# Evolutionary Computation (CEC), 2017 IEEE Congress on. IEEE, 2017.
include("deprecated/CMAES_AEP.jl")

# R.A. Formato Central force optimization:
# a new metaheuristic with applications in applied electromagnetics
# Progress in Electromagnetics Research, PIER 77 (2007), pp. 425-491, 10.2528/PIER07082403
include("deprecated/CFO.jl")

# Brest, J., Maučec, M. S., & Bošković, B. (2017, June).
# Single objective real-parameter optimization: Algorithm jSO.
# In Evolutionary Computation (CEC), 2017 IEEE Congress on (pp. 1311-1318). IEEE.
include("deprecated/jso.jl")


# deprecated functions/optimizers/etc
include("deprecated/deprecated.jl")


end # module
