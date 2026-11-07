MOI.supports(::Optimizer, ::MOI.NLPBlock) = true

function MOI.set(model::Optimizer, ::MOI.NLPBlock, data::Union{Nothing, MOI.NLPBlockData})
    model.nlp_block = data
    return
end

MOI.get(model::Optimizer, ::MOI.NLPBlock) = model.nlp_block

# Nonlinear functions come from two places: the new interface
# (`ScalarNonlinearFunction`, stored in `model.nlp_model`) and the legacy
# `NLPBlock` (`@NLobjective`, `@NLconstraint`). Both evaluators take `x` in MOI
# variable order. Only zero-order evaluations are needed here.
function _nonlinear_sources(model::Optimizer, vars)
    sources = Tuple{Any, Vector{MOI.NLPBoundsPair}}[]
    objective = nothing

    if !isnothing(model.nlp_model)
        backend = MOI.Nonlinear.SparseReverseMode()
        evaluator = MOI.Nonlinear.Evaluator(model.nlp_model, backend, vars)
        MOI.initialize(evaluator, Symbol[])
        bounds = MOI.NLPBlockData(evaluator).constraint_bounds
        isempty(bounds) || push!(sources, (evaluator, bounds))
        isnothing(model.nlp_model.objective) || (objective = evaluator)
    end

    if !isnothing(model.nlp_block)
        block = model.nlp_block
        MOI.initialize(block.evaluator, Symbol[])
        isempty(block.constraint_bounds) ||
            push!(sources, (block.evaluator, block.constraint_bounds))
        # as in other MOI solvers, an NLPBlock objective takes precedence
        block.has_objective && (objective = block.evaluator)
    end

    # FEASIBILITY_SENSE ignores the objective; no objective means f(x) = 0
    model.sense == MOI.FEASIBILITY_SENSE && (objective = nothing)
    return objective, sources
end
