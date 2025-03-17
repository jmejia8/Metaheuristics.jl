include("GA/GA.jl")

# ECA algorithm
include("ECA/ECA.jl")
include("ECA/CECA.jl")

# Differential algorithm
include("DE/DE.jl")

# PSO algorithm
include("PSO/PSO.jl")

# The whale optimization algorithm
# S Mirjalili, A Lewis - Advances in Engineering Software, 2016
include("WOA/WOA.jl")

# Mirjalili, Seyedali, and Amir H. Gandomi.
# "Chaotic gravitational constants for the gravitational search algorithm."
# Applied Soft Computing 53 (2017): 407-419.
include("CGSA/CGSA.jl")

# SA: Simulated Annealing
# Kirkpatrick, S., Gelatt, C.D., & Vecchi, M.P. (1983). Optimization by
# Simulated Annealing. _Science, 220_, 671-680.
include("SA/SA.jl")


# Aritifical Bee colony
include("ABC/ABC.jl")


# genetic algorithm
include("MCCGA/MCCGA.jl")

include("Restart/Restart.jl")
include("SHADE/SHADE.jl")
