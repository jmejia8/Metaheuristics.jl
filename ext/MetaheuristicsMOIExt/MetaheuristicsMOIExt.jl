module MetaheuristicsMOIExt

using Metaheuristics
using SearchSpaces
#using Metaheuristics.LinearAlgebra
import MathOptInterface as MOI

function __init__()
    @static if isdefined(Base, :get_extension)
        setglobal!(Metaheuristics, :Optimizer, Optimizer)
    end
    println("Extension was loaded!")
end

mutable struct Optimizer{T} <: MOI.AbstractOptimizer
    # Problem data.
    variables::MOI.Utilities.VariablesContainer{T}
    starting_values::Vector{Union{Nothing,T}}
    nlp_model::Union{MOI.Nonlinear.Model,Nothing}
    sense::MOI.OptimizationSense

    # Parameters.
    method::Union{Metaheuristics.AbstractAlgorithm,Nothing}
    silent::Bool
    options::Dict{Symbol,Any}

    # Solution attributes.
    results::Union{Nothing,Metaheuristics.State}
    variable_type::Dict
    search_space
    var_dist
    f
    #Union{Nothing,Metaheuristics.MultivariateOptimizationResults}
end

function Optimizer{T}() where {T}
    Optimizer{T}(
                 MOI.Utilities.VariablesContainer{T}(),
                 Union{Nothing,T}[],
                 nothing,
                 MOI.FEASIBILITY_SENSE,
                 nothing,
                 false,
                 Dict{Symbol,Any}(),
                 nothing,
                 Dict(),
                 nothing,
                 nothing,
                 x->x,
                )
end

Optimizer() = Optimizer{Float64}()

MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

function MOI.supports(::Optimizer, ::Union{MOI.ObjectiveSense,MOI.ObjectiveFunction})
    return true
end
MOI.supports(::Optimizer, ::MOI.Silent) = true
function MOI.supports(::Optimizer, p::MOI.RawOptimizerAttribute)
    return p.name == "method" || hasfield(Metaheuristics.Options, Symbol(p.name))
end

function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex})
    # TODO: starting solution
    return false
end


# TODO: add ZeroOne constraints, and integer one
const BOUNDS{T} = Union{MOI.LessThan{T},MOI.GreaterThan{T},MOI.EqualTo{T},MOI.Interval{T}}
const _SETS{T} = Union{MOI.GreaterThan{T},MOI.LessThan{T},MOI.EqualTo{T}}

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.VariableIndex},
    ::Type{<:BOUNDS{T}},
) where {T}
    return true
end


function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{MOI.ZeroOne},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer,
    ::Type{MOI.VariableIndex},
    ::Type{MOI.Integer},
)
    return true
end

function MOI.supports_constraint(
    ::Optimizer{T},
    ::Type{MOI.ScalarNonlinearFunction},
    ::Type{<:_SETS{T}},
) where {T}
    return true
end

MOI.supports_incremental_interface(::Optimizer) = true

function MOI.copy_to(model::Optimizer, src::MOI.ModelLike)
    return MOI.Utilities.default_copy_to(model, src)
end

MOI.get(::Optimizer, ::MOI.SolverName) = "Metaheuristics"

function MOI.set(model::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    model.sense = sense
    return
end
function MOI.set(model::Optimizer, ::MOI.ObjectiveFunction{F}, func::F) where {F}
    # TODO add compat with MOI.VectorNonlinearFunction
    nl = convert(MOI.ScalarNonlinearFunction, func)
    
    if isnothing(model.nlp_model)
        model.nlp_model = MOI.Nonlinear.Model()
    end
    MOI.Nonlinear.set_objective(model.nlp_model, nl)
    nothing
end

function MOI.set(model::Optimizer, ::MOI.Silent, value::Bool)
    model.silent = value
    nothing
end

MOI.get(model::Optimizer, ::MOI.Silent) = model.silent

const TIME_LIMIT = "time_limit"
MOI.supports(::Optimizer, ::MOI.TimeLimitSec) = true
function MOI.set(model::Optimizer, ::MOI.TimeLimitSec, value::Real)
    MOI.set(model, MOI.RawOptimizerAttribute(TIME_LIMIT), Float64(value))
end
function MOI.set(model::Optimizer, attr::MOI.TimeLimitSec, ::Nothing)
    delete!(model.options, Symbol(TIME_LIMIT))
end
function MOI.get(model::Optimizer, ::MOI.TimeLimitSec)
    return get(model.options, Symbol(TIME_LIMIT), nothing)
end

MOI.Utilities.map_indices(::Function, opt::Metaheuristics.AbstractAlgorithm) = opt

function MOI.set(model::Optimizer, p::MOI.RawOptimizerAttribute, value)
    if p.name == "method"
        model.method = value
    else
        model.options[Symbol(p.name)] = value
    end
    return
end

function MOI.get(model::Optimizer, p::MOI.RawOptimizerAttribute)
    if p.name == "method"
        return p.method
    end
    key = Symbol(p.name)
    if haskey(model.options, key)
        return model.options[key]
    end
    error("RawOptimizerAttribute with name $(p.name) is not set.")
end

MOI.get(model::Optimizer, ::MOI.SolveTimeSec) = model.results.overall_time

function MOI.empty!(model::Optimizer)
    MOI.empty!(model.variables)
    empty!(model.starting_values)
    model.nlp_model = nothing
    model.sense = MOI.FEASIBILITY_SENSE
    model.results = nothing
    return
end

function MOI.is_empty(model::Optimizer)
    return MOI.is_empty(model.variables) &&
           isempty(model.starting_values) &&
           isnothing(model.nlp_model) &&
           model.sense == MOI.FEASIBILITY_SENSE
end

function MOI.add_variable(model::Optimizer{T}) where {T}
    push!(model.starting_values, nothing)
    vi = MOI.add_variable(model.variables)
    # Set default type for variables without explicit type constraints
    # This will be overridden if a ZeroOne, Integer, or other constraint is added
    model.variable_type[vi] = :bounds
    return vi
end
function MOI.is_valid(model::Optimizer, index::Union{MOI.VariableIndex,MOI.ConstraintIndex})
    return MOI.is_valid(model.variables, index)
end


function MOI.add_constraint(
    model::Optimizer,
    vi::MOI.VariableIndex,
    set::MOI.ZeroOne,
)
    model.variable_type[vi] = :binary
    # TODO: revise this part to get g(x) <= 0
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    vi::MOI.VariableIndex,
    set::BOUNDS{T},
) where {T}
    model.variable_type[vi] = :bounds
    # TODO: revise this part to get g(x) <= 0
    return MOI.add_constraint(model.variables, vi, set)
end


function MOI.add_constraint(
    model::Optimizer,
    vi::MOI.VariableIndex,
    set::MOI.Integer,
)
    # Mark variable as integer
    model.variable_type[vi] = :int
    return MOI.add_constraint(model.variables, vi, set)
end

function MOI.add_constraint(
    model::Optimizer{T},
    f::MOI.ScalarNonlinearFunction,
    s::_SETS{T},
) where {T}
    if isnothing(model.nlp_model)
        model.nlp_model = MOI.Nonlinear.Model()
    end
    index = MOI.Nonlinear.add_constraint(model.nlp_model, f, s)
    MOI.ConstraintIndex{typeof(f),typeof(s)}(index.value)
end


function starting_value(optimizer::Optimizer{T}, i) where {T}
    if optimizer.starting_values[i] !== nothing
        return optimizer.starting_values[i]
    else
        v = optimizer.variables
        return min(max(zero(T), v.lower[i]), v.upper[i])
    end
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real,Nothing},
)
    MOI.throw_if_not_valid(model, vi)
    model.starting_values[vi.value] = value
    return
end

function _dict_to_search_spaces(dict)

    s = Pair[]
    var_dist = Dict()

    if dict[:binary][:count] > 0
        push!(s, :binary => BitArraySpace(dict[:binary][:count]))
        var_dist[:binary] = dict[:binary][:vars]
    end

    if !isempty(dict[:int][:lb])
        push!(s, :int => BoxConstrainedSpace(dict[:int][:lb], dict[:int][:ub]))
        var_dist[:int] = dict[:int][:vars]
    end

    if !isempty(dict[:bounds][:lb])
        push!(s, :bounds => BoxConstrainedSpace(dict[:bounds][:lb], dict[:bounds][:ub]))
        var_dist[:bounds] = dict[:bounds][:vars]
    end

    if !isempty(dict[:non_rigid][:int][:lb])
        lb = dict[:non_rigid][:int][:lb]
        ub = dict[:non_rigid][:int][:ub]
        push!(s, :int_non_rigid => BoxConstrainedSpace(lb, ub, rigid=false))
        var_dist[:int_non_rigid] = dict[:non_rigid][:int][:vars]
    end

    if !isempty(dict[:non_rigid][:bounds][:lb])
        lb = dict[:non_rigid][:bounds][:lb]
        ub = dict[:non_rigid][:bounds][:ub]
        push!(s, :bounds_non_rigid => BoxConstrainedSpace(lb, ub, rigid=false))
        var_dist[:bounds_non_rigid] = dict[:non_rigid][:bounds][:vars]
    end

    
    MixedSpace(s...), var_dist
end


function _get_search_spaces(model)

    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    lb = model.variables.lower |> copy
    ub = model.variables.upper |> copy

    spaces = Dict(:binary => Dict(:count => 0, :vars => []), 
                  :non_rigid => Dict(
                                     :int => Dict(:lb => Int[], :ub => Int[], :vars => []), 
                                     :bounds => Dict(:lb=>Float64[], :ub=>Float64[], :vars =>[])
                                    ),
                  :int => Dict(:lb => Int[], :ub => Int[], :vars=>[]), 
                  :bounds => Dict(:lb=>Float64[], :ub=>Float64[], :vars=>[])
                 )

    for (a, b, i) in zip(lb, ub, vars)
        # Get variable type, default to :bounds if not set
        k = get(model.variable_type, i, :bounds)
        
        if k === :binary
            spaces[k][:count] += 1
            push!(spaces[k][:vars], i)
            continue
        end

        if isfinite(a) && isfinite(b)
            push!(spaces[k][:lb], a)
            push!(spaces[k][:ub], b)
            push!(spaces[k][:vars], i)
            continue
        end
        
        if isinf(a) && isinf(b)
            # TODO:  study the performance on different initial bound values on
            # unbounded search spaces
            a, b = -10, 10 
        elseif isinf(a)
            a = b - 2max(abs(b), 10)
        else
            b = a + 2max(abs(a), 10)
        end
        push!(spaces[:non_rigid][k][:lb], a)
        push!(spaces[:non_rigid][k][:ub], b)
        push!(spaces[:non_rigid][k][:vars], i)
    end
    spaces
end

function _dict_to_values(solution, search_space, var_dist)
    # TODO improve this part
    x = zeros(Real, Metaheuristics.getdim(search_space))
    for (k, vs) in var_dist
        for (i, v) in enumerate(vs)
            x[v.value] = solution[k][i]
        end
    end
    x
end


function MOI.optimize!(model::Optimizer{T}) where {T}
    #backend = MOI.Nonlinear.ExprGraphOnly() #MOI.Nonlinear.SparseReverseMode()
    backend = MOI.Nonlinear.SparseReverseMode()
    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    evaluator = MOI.Nonlinear.Evaluator(model.nlp_model, backend, vars)
    nlp_data = MOI.NLPBlockData(evaluator)
    objective_scale = model.sense == MOI.MAX_SENSE ? -one(T) : one(T)


    # load parameters
    if isnothing(model.nlp_model)
        error("An objective should be provided to Metaheuristics with `@objective`.")
    end

    nl_constrained = !isempty(nlp_data.constraint_bounds)


    used_features = [:ExprGraph]
    MOI.initialize(evaluator, used_features)


    !model.silent && @info "Load search space"
    search_space_dict = _get_search_spaces(model)
    # get search space and the variable mapping JuMP <--> Metaheuristics
    search_space, var_dist = _dict_to_search_spaces(search_space_dict)
    model.search_space = search_space
    model.var_dist = var_dist

    lb = model.variables.lower 
    ub = model.variables.upper 

    function constraints(x)
        if !nl_constrained
            return #zeros(1), zeros(1)
        end
        
        lc = [b.lower for b in nlp_data.constraint_bounds]
        uc = [b.upper for b in nlp_data.constraint_bounds]


        C = zeros(length(lc));
        MOI.eval_constraint(evaluator, C, x)

        equality_mask =  uc .≈ lc
        inequality_mask =  .!equality_mask
        if count(equality_mask) == 0
            # no equality constraints
            H = zeros(1)
        else
            H = C[equality_mask] - uc[equality_mask]
        end

        # add variable bounds as constraints (due unbounded vars)
        GV = Float64[]
        for (a, b, xv) in zip(lb, ub, x)
            isfinite(a) && push!(GV, a - xv)
            isfinite(b) && push!(GV, xv - b)
        end

        # TODO: The following is hard coded, try to update improve it
        if count(inequality_mask) == 0
            G = zeros(1)
        else
            mask = isfinite.(uc[inequality_mask])
            if any(mask)
                G1 = C[inequality_mask][mask] - uc[inequality_mask][mask] # <= 0
            else
                G1 = Float64[]
            end

            mask = isfinite.(lc[inequality_mask])
            if any(mask)
                G2 = lc[inequality_mask][mask] - C[inequality_mask][mask] # <= 0
            else
                G2 = Float64[]
            end

            
            G = [G1; G2]
        end


        HH = abs.(H) .- 1e-8

        append!(G, GV)
        append!(G, HH)
        

        return G, zeros(1)
    end

    # TODO: add constraints, and multi-objective
    function f(solution)
        x = _dict_to_values(solution, search_space, var_dist)
        fx = objective_scale * MOI.eval_objective(evaluator, x)
        
        C = constraints(x)
        if isnothing(C)
            return fx
        end
        gx, hx = C
        fx, gx, hx
    end


    !model.silent && @info "Load method"

    if isnothing(model.method)
        # TODO add compat with multiobjective optimization
        options = Options(verbose=!model.silent, iterations=10_000,
                          f_calls_limit=Inf)
        algo = Metaheuristics.ECA(adaptive=true)
        model.method = MixedInteger(algo;options)
    end
    model.f = f
    !model.silent && @info "Optimize"

    # TODO initial solutions
    #initial_x = starting_value.(model, eachindex(model.starting_values))
    #options = copy(model.options)
    model.results = optimize(f, search_space, model.method)
    # bounds for constraints


    return
end

function converged(res)
    res.termination_status_code isa Metaheuristics.ConvergenceTermination
end


function MOI.get(model::Optimizer, ::MOI.TerminationStatus)
    if isnothing(model.results)
        return MOI.OPTIMIZE_NOT_CALLED
    end
    
    term = model.results.termination_status_code
    
    if term isa Metaheuristics.ConvergenceTermination
        return MOI.LOCALLY_SOLVED
    elseif term isa Metaheuristics.IterationLimit
        return MOI.ITERATION_LIMIT
    elseif term isa Metaheuristics.EvaluationsLimit
        return MOI.OTHER_LIMIT
    elseif term isa Metaheuristics.TimeLimit
        return MOI.TIME_LIMIT
    elseif term isa Metaheuristics.ObjectiveVarianceLimit
        return MOI.LOCALLY_SOLVED
    elseif term isa Metaheuristics.ObjectiveDifferenceLimit
        return MOI.LOCALLY_SOLVED
    else
        # For other termination reasons (e.g., OtherLimit, UnknownStopReason)
        # Check if we have a feasible solution
        if !isempty(model.results.population)
            best = model.results.best_sol
            if Metaheuristics.is_feasible(best)
                return MOI.LOCALLY_SOLVED
            else
                return MOI.LOCALLY_INFEASIBLE
            end
        end
        return MOI.OTHER_LIMIT
    end
end

function MOI.get(model::Optimizer, ::MOI.RawStatusString)
    sprint(show, model.method)
    #return display(model.results)
end

function MOI.get(model::Optimizer, ::MOI.ResultCount)
    isnothing(model.results) && return 0
    # TODO: check for multi-objective case
    # return length(model.results.population)
    1
end

function MOI.get(model::Optimizer, attr::MOI.PrimalStatus)
    if !(1 <= attr.result_index <= MOI.get(model, MOI.ResultCount()))
        return MOI.NO_SOLUTION
    end
    
    # Check if the best solution is feasible
    if !isempty(model.results.population)
        best = model.results.best_sol
        if Metaheuristics.is_feasible(best)
            return MOI.FEASIBLE_POINT
        else
            return MOI.INFEASIBLE_POINT
        end
    end
    
    return MOI.UNKNOWN_RESULT_STATUS
end
MOI.get(::Optimizer, ::MOI.DualStatus) = MOI.NO_SOLUTION

function MOI.get(model::Optimizer, attr::MOI.ObjectiveValue)
    MOI.check_result_index_bounds(model, attr)
    val = minimum(model.results)
    if model.sense == MOI.MAX_SENSE
        val = -val
    end
    val
end

function MOI.get(model::Optimizer, attr::MOI.VariablePrimal, vi::MOI.VariableIndex)
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, vi)
    solution = Metaheuristics.vec_to_dict(minimizer(model.results), model.search_space)

    x = _dict_to_values(solution, model.search_space, model.var_dist)
    x[vi.value]
end

function MOI.get(
    model::Optimizer{T},
    attr::MOI.ConstraintPrimal,
    ci::MOI.ConstraintIndex{MOI.VariableIndex,<:BOUNDS{T}},
) where {T}
    MOI.check_result_index_bounds(model, attr)
    MOI.throw_if_not_valid(model, ci)
    return Metaheuristics.minimizer(model.results)[ci.value]
end
end # module
