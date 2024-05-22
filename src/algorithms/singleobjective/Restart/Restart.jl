abstract type AbstractRestart <: AbstractParameters end

struct Restart{T} <: AbstractRestart
    alg::T
    every::Int
end

"""
    Restart(optimizer, every=100)

Resets the `optimizer` every specified number of iterations (100 by default).

### Example

```julia-repl
julia> f, bounds, _ = Metaheuristics.TestProblems.rastrigin();

julia> optimize(f, bounds, Restart(ECA(), every=200))
```

### Customization

The restart condition can be updated by overloading the `restart_condition` method:

```julia
function Metaheuristics.restart_condition(status, restart::Restart, information, options)
    st.iteration % params.every == 0
end
```

"""
function Restart(base::Algorithm; every = 100)
    parameters = Restart(base.parameters, every)

    Algorithm(
        parameters,
        information = base.information,
        options = base.options
    )
end

restart_condition(st, params::Restart, args...) = st.iteration % params.every == 0

function restart_population!(st, params::AbstractRestart, problem, info, opts, args...; kargs...)

    if !restart_condition(st, params, info, opts)
        return
    end

    opts.debug && @info "Restarting population..."

    # initialize population (not params.alg)
    st_tmp = State(nothing, [])
    st_new = initialize!(st_tmp, params.alg, problem, info, opts, args...; kargs...)

    # current population is replaced by the new one
    st.population = st_new.population

    opts.debug && @info "Restarted population."
end


function update_state!(st, params::AbstractRestart, problem, info, opts, args...; kargs...)

    # restart population if necessary
    restart_population!(st, params, problem, info, opts, args...; kargs...)

    # perform the optimization step at current population
    update_state!(st, params.alg, problem, info, opts, args...; kargs...)
end


#########################################
# implement each step for AbstractRestart
#########################################
#
function initialize!(st, params::AbstractRestart, problem, info, opts, args...; kargs...)

    initialize!(st, params.alg, problem, info, opts, args...; kargs...)
end


function final_stage!(st, params::AbstractRestart, problem, info, opts, args...; kargs...)

    final_stage!(st, params.alg, problem, info, opts, args...; kargs...)
end


function stop_criteria!(st, params::AbstractRestart, problem, info, opts, args...; kargs...)

    stop_criteria!(st, params.alg, problem, info, opts, args...; kargs...)

    return
end
