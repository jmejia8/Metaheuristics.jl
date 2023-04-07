abstract type AbstractSolution end
abstract type AbstractAlgorithm end
abstract type AbstractParameters end
abstract type AbstractProblem end
abstract type AbstractTermination end

"""
    TerminationStatusCode
An Enum of possible => values for [`State`](@ref).
Possible values:

- `ITERATION_LIMIT`
- `TIME_LIMIT`
- `EVALUATIONS_LIMIT`
- `ACCURACY_LIMIT`
- `OBJECTIVE_VARIANCE_LIMIT`
- `OBJECTIVE_DIFFERENCE_LIMIT`
- `OTHER_LIMIT`
- `UNKNOWN_STOP_REASON`

See also [`termination_status_message`](@ref).
"""
const TerminationStatusCode = AbstractTermination

