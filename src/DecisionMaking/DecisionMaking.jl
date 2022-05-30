abstract type AbstractDecisionMakingMethod end

include("ROI.jl")
include("CompromiseProgramming.jl")

export decisionmaking, dm, best_alternative, ROIArchiving, CompromiseProgramming
export WeightedSum, Tchebysheff, AchievementScalarization

"""
    decisionmaking(fs, w, method)

Perform selected `method` for a given `fs` and weight vector(s) and return the indices
indicating the best alternative(s).
Here, `fs` can be a set of non-dominated solutions (population) or a `State`.
"""
function decisionmaking(population::AbstractArray{<: AbstractMultiObjectiveSolution},
        w,
        method::ROIArchiving
    )

    roiarchiving(population, w, method)
end


decisionmaking(st::State,args...) = decisionmaking(st.population, args...)

const dm = decisionmaking


"""
    best_alternative(population, w, method)

Perform `method` and return the best alternative(s) in `population`.
"""
function best_alternative(
        population::AbstractArray{<: AbstractMultiObjectiveSolution},
        w,
        method::T
    ) where T <: AbstractDecisionMakingMethod

    idx = decisionmaking(population, w, method)
    isempty(idx) && error("Unable finding an alternative using provided method.")

    return population[idx]
end

best_alternative(st::State, args...) = best_alternative(st.population, args...)

