mutable struct Algorithm{T} <: AbstractAlgorithm
    parameters::T
    status::State
    information::Information
    options::Options
    termination::Termination
end

function Algorithm(
        parameters;
        initial_state::State = State(nothing, []),
        information::Information = Information(),
        options::Options = Options(),
        termination = Termination(),
        kargs...
    )
    if !isempty(kargs)
        _pp = collect(kargs)
        p = string.(first.(_pp)) .* " = " .* string.(last.(_pp))
        pp = join(p, ", ")
        @warn "Parameters: `$pp` never used."
    end

    if !(termination isa Termination)
        termination = Termination(checkany = [termination])
    end
    

    Algorithm(parameters, initial_state, information, options, termination)

end

function _print_title(io, title)
    decor = join(repeat('=', length(title)))
    printstyled(io, title, "\n", color=:blue, bold=true)
    printstyled(io, decor, "\n", color=:gray)
end

function Base.show(io::IO, alg::Algorithm)
    _print_title(io, "Algorithm Parameters")
    print(io, "  ")
    Base.show(IOContext(io, :compact => true), alg.parameters)
    println(io, "\n")

    show(io, alg.status)
end

function Base.show(io::IO, parameters::AbstractParameters)
    s = typeof(parameters)

    fns = fieldnames(s)
    all_vals = [getfield(parameters, f) for f in fns]
    # remove those parameters containing empty arrays before showing
    mask = findall(v -> !(v isa Array), all_vals)
    vals = [sprint(show, v) for v in all_vals[mask]]
    str = string(s) * "(" * join(string.(fns[mask]) .* "=" .* vals, ", ") * ")"

    print(io, str)
end

termination_status_message(alg::Algorithm) = termination_status_message(alg.status)
should_stop(algorithm::AbstractAlgorithm) = algorithm.status.stop
function get_result(algorithm::AbstractAlgorithm)
    if isnothing(algorithm.status.best_sol)
        error("First optimize a function: `optimize!(f, bounds, method)`")
    end
    
    algorithm.status
end



iscompatible(search_space, algorithm) = false
iscompatible(search_space::BoxConstrainedSpace, algorithm::AbstractParameters) = true
iscompatible(search_space::AbstractMatrix, algorithm::AbstractParameters) = true

function check_compat(search_space, algorithm::AbstractParameters)
    if iscompatible(search_space, algorithm)
        return
    end
    
    a = typeof(algorithm)
    s = typeof(search_space)
    error("Implemented metaheuristic `$a` is not compatible with search space `$s`. Try using bounds.")
end

