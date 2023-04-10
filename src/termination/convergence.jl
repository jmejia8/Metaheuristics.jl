abstract type ConvergenceTermination <: AbstractTermination end

function termination_status_message(criterion::ConvergenceTermination)
    "Due to Convergence Termination criterion."
end

function stop_check(status::State, criterion::ConvergenceTermination; report=true)
    population = status.population
    # nothing to do for empty populations
    isempty(population) && return false
    # check multi-objective problem
    fval(first(status.population)) isa Array && return false

    feasible = is_feasible.(population)

    # check whether feasibility is ratio over 30%
    count(feasible) / length(population) < 0.3 && return false

    stop = stop_check(population[feasible], criterion)
    
    report && stop && (status.termination_status_code = criterion)
    
    stop
end

"""
    RelativeFunctionConvergence(;ftol)

Check for a small relative difference between the function values of the
current population with the largest and smallest function values.
"""
Base.@kwdef struct RelativeFunctionConvergence <: ConvergenceTermination
    ftol::Float64 = eps()
end


function stop_check(population::AbstractVector, criterion::RelativeFunctionConvergence)

    m, M = extrema(fvals(population))
    abs(M - m) / max(abs(M), eps()) <= criterion.ftol
end

"""
    SmallStandardDeviation(;ftol)

Termination requires a small standard deviation of the function values of
the population.
"""
Base.@kwdef struct SmallStandardDeviation <: ConvergenceTermination
    ftol::Float64 = 1e-6
end

function stop_check(population::AbstractVector, criterion::SmallStandardDeviation)
    fx = fvals(population)
    std(fx) <= criterion.ftol
end

"""
    AbsoluteFunctionConvergence(;ftol)

Termination requires a small absolute difference between the function
values of the population with the largest and smallest function values.
"""
Base.@kwdef struct AbsoluteFunctionConvergence <: ConvergenceTermination
    ftol::Float64 = 0.0
end

function stop_check(population::AbstractVector, criterion::AbsoluteFunctionConvergence)
    m, M = extrema(fvals(population))
    abs(M - m) <= criterion.ftol
end

"""
    RelativeParameterConvergence(;xtol)

Termination requires a small relative parameter difference between the
vertices with the largest and smallest function values.
"""
Base.@kwdef struct RelativeParameterConvergence <: ConvergenceTermination
    xtol::Float64 = 1e-8
end

function stop_check(population::AbstractVector, criterion::RelativeParameterConvergence)
    
    fx = fvals(population)
    xhi = population[argmax(fx)] |> get_position
    xlo = population[argmin(fx)] |> get_position
    # check for numerical representation
    if eltype(xhi) <: Number
        a = max(abs.(xlo), abs.(xhi))
        return maximum(abs.(xhi -xlo))/max(maximum(a), eps()) <= criterion.xtol
    end
    xhi == xlo
end

struct CheckConvergence <: ConvergenceTermination
    criteria::Vector
end

function CheckConvergence(;
        f_tol_rel = eps(),
        f_tol_abs = 0.0,
        x_tol = 1e-8,
        criteria = [
                    AbsoluteFunctionConvergence(f_tol_abs),
                    RelativeFunctionConvergence(f_tol_rel),
                    SmallStandardDeviation(),
                    RelativeParameterConvergence(x_tol),
                   ]
    )
    CheckConvergence(criteria)
end

function stop_check(population::AbstractVector, criteria::CheckConvergence)
    all(stop_check(population, criterion) for criterion in criteria.criteria)
end
