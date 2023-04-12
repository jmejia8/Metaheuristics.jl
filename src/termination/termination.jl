include("stop_status_codes.jl")
include("convergence.jl")
include("budget.jl")
include("default.jl")

Base.@kwdef struct Termination <: AbstractTermination
    checkall::Vector = Any[]
    checkany::Vector = Any[]
end

function stop_check(status, termination::Termination)
    checkall = termination.checkall
    # all criteria are satisfied
    !isempty(checkall) && return all(stop_check(status, c) for c in checkall)

    # check if any other does
    any(stop_check(status, c) for c in termination.checkany)
end


