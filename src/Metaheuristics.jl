__precompile__()
module Metaheuristics

# ECA algorithm
export eca, diffEvolution, pso, WOA, GOA, GSA, CGSA, SA, CMAES_AEP
export CFO, jso, DE, ABC, PSO

# v2 items
export optimize, ECA, Options, State, Information, Engine, Problem

include("externals.jl")
include("stop.jl")
include("structures.jl")
include("operators.jl")
include("display.jl")
include("algorithm.jl")
include("tools.jl")

include("optimize.jl")

# ECA algorithm
include("eca.jl")

# Differential algorithm
include("diffEvolution.jl")

# PSO algorithm
include("pso.jl")

# The whale optimization algorithm
# S Mirjalili, A Lewis - Advances in Engineering Software, 2016
include("WOA.jl")

# Grasshopper optimisation algorithm: Theory and application
# S Saremi, S Mirjalili, A Lewis - Advances in Engineering Software, 2017
include("GOA.jl")

# GSA: a gravitational search algorithm
# E Rashedi, H Nezamabadi-Pour, S Saryazdi - Information sciences, 2009
include("GSA.jl")

# Mirjalili, Seyedali, and Amir H. Gandomi. 
# "Chaotic gravitational constants for the gravitational search algorithm." 
# Applied Soft Computing 53 (2017): 407-419.
include("CGSA.jl")

# SA: Simulated Annealing
# Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
# Simulated Annealing. _Science, 220_, 671-680.
include("SA.jl")

# Li, Zhenhua, and Qingfu Zhang. 
# "An efficient rank-1 update for Cholesky CMA-ES using auxiliary evolution path." 
# Evolutionary Computation (CEC), 2017 IEEE Congress on. IEEE, 2017.
include("CMAES_AEP.jl")

# R.A. Formato Central force optimization:
# a new metaheuristic with applications in applied electromagnetics
# Progress in Electromagnetics Research, PIER 77 (2007), pp. 425-491, 10.2528/PIER07082403
include("CFO.jl")

# Brest, J., Maučec, M. S., & Bošković, B. (2017, June).
# Single objective real-parameter optimization: Algorithm jSO.
# In Evolutionary Computation (CEC), 2017 IEEE Congress on (pp. 1311-1318). IEEE.
include("jso.jl")


include("ABC.jl")

include("deprecated.jl")


end # module
