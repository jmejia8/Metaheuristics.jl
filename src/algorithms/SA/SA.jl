# Based on MATLAB code of Héctor Corte
# B.Sc. in physics 2010


include("new_solution.jl")

mutable struct SA <: AbstractParameters
    N::Int
    x_initial::Vector{Float64}
    tol_fun::Float64
    x::Vector{Float64}
    fx::Float64
end

"""

```julia
    SA(;
        x_initial::Vector = zeros(0),
        N::Int = 500,
        tol_fun::Real= 1e-4,
        information = Information(),
        options = Options()
    )
``` 

Parameters for the method of Simulated Annealing (Kirkpatrick et al., 1983).

Parameters:

- x_intial: Inital solution. If empty, then SA will generate a random one within the bounds.
- N: The number of test points per iteration.
- tol_fun: tolerance value for the Metropolis condition to accept or reject the test point as current point.
 

# Example


```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], SA())
+=========== RESULT ==========+
  iteration: 60
    minimum: 5.0787e-68
  minimizer: [-2.2522059499734615e-34, 3.816133503985569e-36, 6.934348004465088e-36]
    f calls: 29002
 total time: 0.0943 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], SA(N = 100, x_initial = [1, 0.5, -1]))
+=========== RESULT ==========+
  iteration: 300
    minimum: 1.99651e-69
  minimizer: [4.4638292404181215e-35, -1.738939846089388e-36, -9.542441152683457e-37]
    f calls: 29802
 total time: 0.0965 s
+============================+
```


"""
function SA(;
        x_initial::Vector = zeros(0),
        N::Int = 500,
        tol_fun::Real= 1e-4,
        kargs...
    )
    
    parameters = SA(N, x_initial, tol_fun, zeros(0), Inf)

    Algorithm( parameters; kargs...)

    
end

function initialize!(
    status,
    parameters::SA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
   )


    l = problem.search_space.lb
    u = problem.search_space.ub
    rng = options.rng

    if isempty(parameters.x_initial)
        parameters.x_initial = l .+ (u .- l) .* rand(rng, length(u))
	end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000length(u)
    end

    if options.iterations == 0
        options.iterations = options.f_calls_limit ÷ parameters.N
    end

	# the current point and fx=f(x)
	x = parameters.x_initial
    best_sol = create_solution(x, problem)
	fx = best_sol.f
    status = State(best_sol, [best_sol])
    parameters.x = x
    parameters.fx = fx
    status.f_calls = 1

    return status

end


function update_state!(
    status::State{xf_indiv},
    parameters::SA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )

    nevals = status.f_calls
    max_evals = options.f_calls_limit
    N = parameters.N
    TolFun = parameters.tol_fun
    rng = options.rng

    # T is the inverse of temperature.
    T = nevals / max_evals 
    μ = 10.0 ^( 100T )    

    l = problem.search_space.lb
    u = problem.search_space.ub

    # For each temperature we take 500 test points to simulate reach termal
    # equilibrium.
    for i = 1:N        
        # We generate new test point using newSol function      
        dx = newSol(2rand(rng, length(parameters.x)) .- 1.0 , μ) .* (u-l)

        # the test point and fx1=f(x1)
        x1 = parameters.x + dx

        # Next step is to keep solution within bounds
        #x1 = (x1 .< l).*l+(l .<= x1).*(x1 .<= u).*x1+(u .< x1).*u			
        reset_to_violated_bounds!(x1, problem.search_space)
        fx1 = evaluate(x1, problem)

        status.f_calls += 1


        df  = fx1 - parameters.fx

        # If the function variation,df, is <0 we take test point as current
        # point. And if df>0 we use Metropolis condition to accept or
        # reject the test point as current point.
        if (df < 0 || rand(rng) < exp(-T*df/(abs(parameters.fx)) / TolFun))
            parameters.x = x1
            parameters.fx= fx1
        end        

        # If the current point is better than current solution, we take
        # current point as cuyrrent solution.       
        if fx1 < status.best_sol.f
            status.best_sol.x = x1
            status.best_sol.f = fx1
        end

        stop_criteria!(status, parameters, problem, information, options)

        if status.stop
            break
        end
    end


end


function final_stage!(
    status::State,
    parameters::SA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
    )
    status.final_time = time()
    status.population = [status.best_sol]
end

