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
    codes = Dict(
                 ITERATION_LIMIT => "Maximum number of iterations exceeded.",
                 TIME_LIMIT => "Maximum time exceeded.",
                 EVALUATIONS_LIMIT => "Maximum objective function calls exceeded.",
                 ACCURACY_LIMIT => "The desired accuracy was obtained.",
                 OBJECTIVE_VARIANCE_LIMIT => "Small variance of the objective function.",
                 OBJECTIVE_DIFFERENCE_LIMIT =>"Small difference of objective function values.",
                 OTHER_LIMIT => "Other stopping criteria.",
                 UNKNOWN_STOP_REASON => "Unknown stop reason."
    )
    try
        return codes[status_code]
    catch
        throw("Illegal value $status_code")
    end
end

