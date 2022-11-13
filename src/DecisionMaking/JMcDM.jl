const mcdm = JMcDM.mcdm
const MCDMSetting = JMcDM.MCDMSetting
const summary = JMcDM.summary


export mcdm
export MCDMSetting
export summary

export best_alternative


"""
    mcdm(fs, w, method)

Perform selected `method` for a given `fs` and weight vector.
Here, `fs` can be a set of non-dominated solutions (population), a `State`
or a decision matrix.

Also, `method` can be selected from `JMcDM` package.

Supported McDM methods:

* `ArasMethod`
* `CocosoMethod`
* `CodasMethod`
* `CoprasMethod`
* `EdasMethod`
* `ElectreMethod`
* `GreyMethod`
* `MabacMethod`
* `MaircaMethod`
* `MooraMethod`
* `SawMethod`
* `TopsisMethod` **(default method)**
* `VikorMethod`
* `WPMMethod`
* `WaspasMethod`
* `MarcosMethod`

### Example 1:

Performing McDM using a population.

```julia-repl
julia> using JMcDM

julia> _, _, population = Metaheuristics.TestProblems.ZDT1();

julia> dm = mcdm(population, [0.5, 0.5], TopsisMethod());

julia> population[dm.bestIndex]
(f = [0.5353535353535354, 0.2683214262030523], g = [0.0], h = [0.0], x = [5.354e-01, 0.000e+00, …, 0.000e+00])
```

### Example 2:

Performing McDM using results from metaheuristic.

```julia-repl
julia> using JMcDM

julia> f, bounds, _ = Metaheuristics.TestProblems.ZDT1();

julia> res = optimize(f, bounds, NSGA2());

julia> dm = mcdm(res, [0.5, 0.5], TopsisMethod());

julia> res.population[dm.bestIndex]
(f = [0.32301132058506055, 0.43208538139854685], g = [0.0], h = [0.0], x = [3.230e-01, 1.919e-04, …, 1.353e-04])
```
"""
function JMcDM.mcdm(
    f::AbstractMatrix{<:AbstractFloat}, # objective functions by col
    w,
    args...
    )

    fns = [minimum for i in 1:size(f, 1)]
    JMcDM.mcdm(JMcDM.DataFrame(f, :auto), w, fns, args...)
end

function JMcDM.mcdm(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    fs = fvals(population) # Pareto solutions
    mcdm(fs, args...)
end

function JMcDM.mcdm(st::State,args...)
    mcdm(st.population, args...)
end

function JMcDM.MCDMSetting(f::AbstractMatrix{<: AbstractFloat}, weights)
    fns = [minimum for i in 1:size(f, 1)]
    JMcDM.MCDMSetting(JMcDM.DataFrame(f, :auto), weights, fns)
end


function JMcDM.MCDMSetting(status::State,args...)
    MCDMSetting(status.population, args...)
end

function JMcDM.MCDMSetting(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    MCDMSetting(fvals(population), args...)
end

function JMcDM.summary(
    f::AbstractMatrix{<: AbstractFloat},
    w::Array{Float64,1}, 
    methods::Array{Symbol,1}
    )

    fns = [minimum for i in 1:size(f, 1)]
    JMcDM.summary(JMcDM.DataFrame(f, :auto), w, fns, methods)
end


function JMcDM.summary(population::AbstractArray{<: AbstractMultiObjectiveSolution},args...)
    summary(fvals(population), args...)
end


function JMcDM.summary(st::State,args...)
    summary(st.population, args...)
end


function decisionmaking(
    f::AbstractMatrix{<:AbstractFloat}, # objective functions by col
    w::AbstractVector{<:AbstractFloat},
    method::T
    ) where T <: JMcDM.MCDMMethod

    result = mcdm(f, w, method)

    if result.bestIndex isa Tuple
        return [result.bestIndex...]
    end

    result.bestIndex
end

function decisionmaking(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    decisionmaking(fvals(population), args...)
end


function decisionmaking(st::State, args...)
    decisionmaking(st.population, args...)
end


"""
    best_alternative(res, w, method)

Perform McDM using results from metaheuristic and return best alternative in `res.population`.

### Example

```julia-repl
julia> f, bounds, _ = Metaheuristics.TestProblems.ZDT1();

julia> res = optimize(f, bounds, NSGA2());

julia> best_sol = best_alternative(res, [0.5, 0.5], TopsisMethod())
(f = [0.32301132058506055, 0.43208538139854685], g = [0.0], h = [0.0], x = [3.230e-01, 1.919e-04, …, 1.353e-04])
```
"""
function best_alternative(
    population::AbstractArray{<: AbstractMultiObjectiveSolution},
    w::AbstractVector{<:AbstractFloat},
    method::T
    ) where T <: JMcDM.MCDMMethod

    idx = decisionmaking(population, w, method)
    isempty(idx) && error("Unable finding an alternative using provided method.")
 
    population[idx]
end

function best_alternative(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    best_alternative(fvals(population), args...)
end


function best_alternative(st::State, args...)
    best_alternative(st.population, args...)
end

best_alternative(status::State, args...) = best_alternative(status.population,args...)




