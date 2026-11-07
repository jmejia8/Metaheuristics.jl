# Start values may be partial: see `_start_solutions`.
function MOI.supports(::Optimizer, ::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex})
    return true
end

function MOI.set(
    model::Optimizer,
    ::MOI.VariablePrimalStart,
    vi::MOI.VariableIndex,
    value::Union{Real, Nothing},
)
    MOI.throw_if_not_valid(model, vi)
    model.starting_values[vi.value] = value
    return
end

function MOI.get(model::Optimizer, ::MOI.VariablePrimalStart, vi::MOI.VariableIndex)
    MOI.throw_if_not_valid(model, vi)
    return model.starting_values[vi.value]
end

# Population size of `method`, or `nothing` when it is chosen in `initialize!`
# (e.g. ECA, DE and PSO with the default N = 0).
function _population_size(method)
    p = _base_parameters(method.parameters)
    N = hasproperty(p, :N) ? p.N : 0
    return N isa Integer && N > 0 ? Int(N) : nothing
end

function _n_start_seeds(model::Optimizer, method, D, partial::Bool)
    partial || return 1 # every value is given: more copies would be duplicates
    N = _population_size(method)
    k = something(model.start_seeds, ceil(Int, START_SEEDS_FRACTION * something(N, 10D)))
    # keep the population size (MOEA/D and NSGA-III tie it to their weights)
    return isnothing(N) ? k : clamp(k, 1, N)
end

_is_given(v) = !isnothing(v) && isfinite(v)

# `k` positions in the MixedSpace order used by `MixedInteger`, built from the start
# values. Given values are kept (rounded for integer and binary variables, clamped
# to rigid bounds), and missing ones are sampled uniformly in the search space. A
# partial start thus defines a subspace, and each position is a different point of
# it. Also returns how many given values had to be rounded or clamped.
function _start_positions(model::Optimizer, search_space, var_dist, k, rng)
    X = [Float64[] for _ in 1:k]
    n_changed = 0
    for key in search_space.key_order
        space = search_space.domain[key]
        binary = space isa BitArraySpace
        discrete = binary || eltype(space.lb) <: Integer
        for (j, vi) in enumerate(var_dist[key])
            a, b = binary ? (0.0, 1.0) : (Float64(space.lb[j]), Float64(space.ub[j]))
            v = model.starting_values[vi.value]
            if _is_given(v)
                w = discrete ? round(Float64(v)) : Float64(v)
                (binary || space.rigid) && (w = clamp(w, a, b))
                n_changed += w != v
                foreach(x -> push!(x, w), X)
            elseif discrete
                foreach(x -> push!(x, Float64(rand(rng, Int(a):Int(b)))), X)
            else
                foreach(x -> push!(x, a + (b - a) * rand(rng)), X)
            end
        end
    end
    return X, n_changed
end

# Add the solutions built from the start values (if any) to the initial state of
# `method`; the base algorithm completes the population with random solutions.
function _set_start_solutions!(method, model::Optimizer, f, search_space, var_dist)
    starts = model.starting_values
    any(_is_given, starts) || return
    partial = !all(_is_given, starts)
    D = Metaheuristics.getdim(search_space)
    k = _n_start_seeds(model, method, D, partial)
    rng = method.options.rng
    X, n_changed = _start_positions(model, search_space, var_dist, k, rng)
    if n_changed > 0 && !model.silent
        @warn "$n_changed start value(s) rounded or clamped to the variable type and bounds."
    end
    for x in X
        # `vec_to_dict` also writes the rounded values back into `x`
        fx = f(Metaheuristics.vec_to_dict(x, search_space))
        sol = Metaheuristics.create_child(x, fx; ε = method.options.h_tol)
        Metaheuristics.set_user_solutions!(method, sol)
    end
    return
end
