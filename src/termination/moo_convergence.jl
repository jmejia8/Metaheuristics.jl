abstract type MOOConvergenceTermination <: AbstractTermination end


Base.@kwdef mutable struct RobustConvergence <: ConvergenceTermination
    ftol::Float64 = 1e-3
    period::Int = 10
    _last_population = []
    _iter::Int = 0
end

function _check_extremas(fa, fb, tol)
    maximum(abs.(fb - fa)) <= tol
end

function stop_check(f_current, f_last, criterion::RobustConvergence)
    tol = criterion.ftol
    # compute extrema
    fmin = ideal(f_current)
    fmax = nadir(f_current)
    # normalization denominator
    Δ = fmax - fmin
    Δ[Δ .≈ 0] .= one(eltype(Δ))
    
    # if extrema are not in tolerance, then don't stop
    !_check_extremas(fmin ./ Δ, ideal(f_last) ./ Δ, tol) && return false
    !_check_extremas(fmax ./ Δ, nadir(f_last) ./ Δ, tol) && return false 
    
    
    
    # normalize fronts
    f_current = (f_current .- fmin') ./ Δ'
    f_last = (f_last .- fmin') ./ Δ'
    # stop if extrema and now IGD is under tolerance
    indicator = PerformanceIndicators.igd(f_last, f_current)
    indicator <= tol
end

function stop_check(population::AbstractVector, criterion::RobustConvergence)
    last_population = criterion._last_population
    if isempty(last_population)
        criterion._last_population = deepcopy(population)
    end

    criterion._iter += 1
    criterion._iter % criterion.period != 0 && return false

    stop = stop_check(fvals(population), fvals(last_population), criterion)
    criterion._last_population = deepcopy(population)
    stop
end

