# Copy the user options (`RawOptimizerAttribute`, `TimeLimitSec`, `Silent`)
# into the options of the method used for solving.
function _apply_options!(opts::Metaheuristics.Options, model::Optimizer)
    for (key, value) in model.options
        setfield!(opts, key, value)
    end
    # a new seed needs a new generator, unless the user also gave one
    if haskey(model.options, :seed) && !haskey(model.options, :rng)
        opts.rng = Metaheuristics.default_rng_mh(opts.seed)
    end
    # Silent = false keeps the verbosity chosen by the user
    model.silent && (opts.verbose = false)
    return opts
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
        push!(s, :int_non_rigid => BoxConstrainedSpace(lb, ub, rigid = false))
        var_dist[:int_non_rigid] = dict[:non_rigid][:int][:vars]
    end

    if !isempty(dict[:non_rigid][:bounds][:lb])
        lb = dict[:non_rigid][:bounds][:lb]
        ub = dict[:non_rigid][:bounds][:ub]
        push!(s, :bounds_non_rigid => BoxConstrainedSpace(lb, ub, rigid = false))
        var_dist[:bounds_non_rigid] = dict[:non_rigid][:bounds][:vars]
    end

    MixedSpace(s...), var_dist
end

# Return a message if the variable bounds admit no point, otherwise `nothing`.
function _check_bounds(model)
    lb = model.variables.lower
    ub = model.variables.upper
    for vi in MOI.get(model.variables, MOI.ListOfVariableIndices())
        a, b = lb[vi.value], ub[vi.value]
        k = get(model.variable_type, vi, :bounds)
        if k === :binary
            a, b = max(a, 0), min(b, 1)
        end
        if k !== :bounds
            a, b = ceil(a), floor(b)
        end
        if a > b
            kind = k === :bounds ? "" : k === :int ? "integer " : "binary "
            return "Infeasible bounds: $(kind)variable $(vi.value) has no value " *
                "in [$(lb[vi.value]), $(ub[vi.value])]."
        end
    end
    return nothing
end

function _get_search_spaces(model)

    vars = MOI.get(model.variables, MOI.ListOfVariableIndices())
    lb = model.variables.lower |> copy
    ub = model.variables.upper |> copy

    spaces = Dict(
        :binary => Dict(:count => 0, :vars => []),
        :non_rigid => Dict(
            :int => Dict(:lb => Int[], :ub => Int[], :vars => []),
            :bounds => Dict(:lb => Float64[], :ub => Float64[], :vars => []),
        ),
        :int => Dict(:lb => Int[], :ub => Int[], :vars => []),
        :bounds => Dict(:lb => Float64[], :ub => Float64[], :vars => []),
    )

    for (a, b, i) in zip(lb, ub, vars)
        # Get variable type, default to :bounds if not set
        k = get(model.variable_type, i, :bounds)

        if k === :binary
            spaces[k][:count] += 1
            push!(spaces[k][:vars], i)
            continue
        end

        if k === :int
            # Too large bounds cannot be stored as Int: handle them as infinite
            # here (they are still enforced as penalty constraints).
            abs(a) < INT_BOUND_LIMIT || (a = -Inf)
            abs(b) < INT_BOUND_LIMIT || (b = Inf)
        end

        if isfinite(a) && isfinite(b)
            if k === :int
                # fractional bounds: tightest integer box
                a, b = ceil(Int, a), floor(Int, b)
            end
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
        if k === :int
            a, b = ceil(Int, a), floor(Int, b)
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
