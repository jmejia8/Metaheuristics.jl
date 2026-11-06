# Equality constraints h(x) = 0 are handled as |h(x)| - tol <= 0, where tol is
# the `h_tol` of the method options, or this value when `h_tol` is zero.
const DEFAULT_H_TOL = 1.0e-6

const START_SEEDS = "start_seeds"

# TODO: add ZeroOne constraints, and integer one
const BOUNDS{T} = Union{
    MOI.LessThan{T},
    MOI.GreaterThan{T},
    MOI.EqualTo{T},
    MOI.Interval{T},
}
const _SETS{T} = Union{MOI.GreaterThan{T}, MOI.LessThan{T}, MOI.EqualTo{T}}

const TIME_LIMIT = "time_limit"

# Default number of initial solutions built from a partial start, as a fraction of
# the population size N (N is estimated as 10 D when the method chooses it later).
const START_SEEDS_FRACTION = 0.1

# Integer bounds at or beyond this value are treated as infinite in the search
# space (an `Int` box cannot hold them, and floats lose integer precision).
const INT_BOUND_LIMIT = 2.0^52

# Algorithms that need a vector objective.
const _MO_PARAMETERS = Union{
    Metaheuristics.NSGA2,
    Metaheuristics.NSGA3,
    Metaheuristics.SMS_EMOA,
    Metaheuristics.SPEA2,
    Metaheuristics.MOEAD_DE,
    Metaheuristics.CCMO,
}

# Stops meaning "the search settled" (not "the budget ran out").
const _SETTLED = Union{
    Metaheuristics.ConvergenceTermination,
    Metaheuristics.ObjectiveVarianceLimit,
    Metaheuristics.ObjectiveDifferenceLimit,
    Metaheuristics.AccuracyLimit,
}
