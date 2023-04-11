abstract type MOOConvergenceTermination <: AbstractTermination end


Base.@kwdef mutable struct RobustConvergence <: ConvergenceTermination
    ftol::Float64 = 1e-3
    period::Int = 10
    _last_population = []
    _iter::Int = 0
end

function check_extremas(fa, fb, tol)
    @show maximum(abs.(fb - fa))
    maximum(abs.(fb - fa)) < tol
end

function stop_check(f_current, f_last, criterion::RobustConvergence)
    tol = criterion.ftol
    # compute extrema
    fmin = ideal(f_current)
    fmax = nadir(f_current)
    # normalization denominator
    Δ = fmax - fmin
    Δ[dnom .≈ 0] .= 1
    
    # if extrema are not in tolerance, then don't stop
    !check_extremas(fmin ./ Δ, ideal(f_last) ./ Δ, tol) && return false
    !check_extremas(fmax ./ Δ, nadir(f_last) ./ Δ, tol) && return false 
    
    # normalize fronts
    f_current = (fs - fmin) ./ Δ
    f_last = (fvals(last_population) - fmin) ./ Δ
    # stop if extrema and now IGD is under tolerance
    indicator = PerformanceIndicators.igd(f_prev, f_current)
    @show indicator
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

