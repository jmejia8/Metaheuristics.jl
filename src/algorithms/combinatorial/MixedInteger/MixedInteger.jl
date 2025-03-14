# structure with algorithm parameters
mutable struct MixedInteger{T} <: AbstractParameters
    base::T
    problem
end

include("utils.jl")

_num_to_vec(v::AbstractDict) = v
iscompatible(search_space::MixedSpace, algorithm::MixedInteger) = true

"""
    MixedInteger(algorithm; information = Information(), options = Options())

Create a wrapper to enable an `algorithm` to handle mixed-integer problems.
When `MixedInteger` is used, the objective function will receive a dictionary with values
in corresponding search space.

# Example

```julia

# objective function
f(solution) = sum(solution[:integer]) * (1 .- sum(solution[:continuous])^2 )

# search spaces
integer_space = BoxConstrainedSpace(-10ones(Int, 6), 10ones(Int, 6))
continuous_space = BoxConstrainedSpace(-ones(3), ones(3))
# mix spaces to form a mixed space
mixed_space = MixedSpace(:integer => integer_space, :continuous => continuous_space)

# wrap ECA to handle mixed-integer variables
mip_eca = MixedInteger(ECA(N = 100), options = Options(seed = 1))

result = optimize(f, mixed_space, mip_eca)

# convert resulting minimizer (which is a vector) into a dictionary
solution = Metaheuristics.vec_to_dict(minimizer(result), mixed_space)
display(solution)
```

**output**

```julia
Dict{Symbol, Vector} with 2 entries:
  :continuous => [-1.0, -1.0, -1.0]
  :integer    => [10, 10, 10, 10, 10, 10]
```

"""
function MixedInteger(base; information = Information(), options = Options())
    parameters = MixedInteger(base.parameters, nothing)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )
end


function initialize!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    # introduce a helper objective function that translates vectors into dictionaries
    f(v) = begin
        d = vec_to_dict(v, problem.search_space)
        problem.f(d)
    end

    # create bounds from mixed
    bounds = _mixed_to_continous_space(problem.search_space)
    parameters.problem = Problem(f, bounds)


    status = initialize!(status, parameters.base, parameters.problem, information, options, args...; kargs...)

    #st = gen_initial_state(parameters.problem,parameters.base,information,options,status)

    problem.f_calls = parameters.problem.f_calls

    status
end


function update_state!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )

    update_state!(status, parameters.base, parameters.problem, information, options, args...; kargs...)
    problem.f_calls = parameters.problem.f_calls
end

function final_stage!(
        status,
        parameters::MixedInteger,
        problem,
        information,
        options,
        args...;
        kargs...
    )
    
    # TODO: convert vectors in population into dictionaries
end


