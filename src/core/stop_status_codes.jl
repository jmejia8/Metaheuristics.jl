"""
    TerminationStatusCode
An Enum of possible return values for [`State`](@ref).
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
@enum(TerminationStatusCode,
      ITERATION_LIMIT, 
      TIME_LIMIT, 
      EVALUATIONS_LIMIT, 
      ACCURACY_LIMIT, 
      OBJECTIVE_VARIANCE_LIMIT,
      OBJECTIVE_DIFFERENCE_LIMIT,
      OTHER_LIMIT,
      UNKNOWN_STOP_REASON,
     )

"""
    termination_status_message(status)

Return a string of the message related to the `status`.

See [`TerminationStatusCode`](@ref)

### Example:

```julia-repl
julia> termination_status_message(Metaheuristics.ITERATION_LIMIT)
"Maximum number of iterations exceeded."

julia> termination_status_message(optimize(f, bounds))
"Maximum number of iterations exceeded."

julia> termination_status_message(ECA())
"Unknown stop reason."
```
"""
function termination_status_message(status_code::TerminationStatusCode)
    if status_code == ITERATION_LIMIT
        return "Maximum number of iterations exceeded."
    elseif status_code == TIME_LIMIT
        return "Maximum time exceeded."
    elseif status_code == EVALUATIONS_LIMIT
        return "Maximum objective function calls exceeded."
    elseif status_code == ACCURACY_LIMIT
        return "The desired accuracy was obtained."
    elseif status_code == OBJECTIVE_VARIANCE_LIMIT
        return "Small variance of the objective function."
    elseif status_code == OBJECTIVE_DIFFERENCE_LIMIT
        return "Small Difference of the max and min value of the objective."
    elseif status_code == OTHER_LIMIT
        return "Other stopping criteria."
    elseif status_code == UNKNOWN_STOP_REASON
        return "Unknown stop reason."
    end

    throw("Illegal value $status_code")
end

