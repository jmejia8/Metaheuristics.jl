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
| Iter.: 60
| f(x) = 2.84574e-73
| solution.x = [-5.307880224731971e-37, -5.183298967486749e-38, 1.2301984439451926e-38]
| f calls: 29502
| Total time: 0.0465 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], SA(N = 100, x_initial = [1, 0.5, -1]))
+=========== RESULT ==========+
| Iter.: 300
| f(x) = 1.29349e-70
| solution.x = [-7.62307964668667e-36, 8.432089040013441e-36, -3.7077496015659554e-37]
| f calls: 29902
| Total time: 0.0466 s
+============================+
```


"""
function SA(;
        x_initial::Vector = zeros(0),
        N::Int = 500,
        tol_fun::Real= 1e-4,
        information = Information(),
        options = Options()
    )
    
    parameters = SA(N, x_initial, tol_fun, zeros(0), Inf)

    Algorithm(
        parameters,
        information = information,
        options = options,
    )

    
end

function initialize!(
    parameters::SA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
   )


    l = view( problem.bounds, 1,:)
    u = view( problem.bounds, 2,:)

    if isempty(parameters.x_initial)
        parameters.x_initial = l .+ (u .- l) .* rand(length(u))
	end

    if options.f_calls_limit == 0
        options.f_calls_limit = 10000length(u)
    end

    if options.iterations == 0
        options.iterations = options.f_calls_limit ÷ parameters.N
    end

	# the current point and fx=f(x)
	x = parameters.x_initial
	fx= problem.f(x)
    best_sol = generateChild(x, fx)
    status = State(best_sol, [best_sol])
    parameters.x = x
    parameters.fx = fx
    status.f_calls = 1

    return status

end


function update_state!(
    status::State,
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

    # T is the inverse of temperature.
    T = nevals / max_evals 
    μ = 10.0 ^( 100T )    

    l = view( problem.bounds, 1,:)
    u = view( problem.bounds, 2,:)

    # For each temperature we take 500 test points to simulate reach termal
    # equilibrium.
    for i = 1:N        
        # We generate new test point using newSol function      
        dx = newSol(2rand(length(parameters.x)) .- 1.0 , μ) .* (u-l)

        # the test point and fx1=f(x1)
        x1 = parameters.x + dx

        # Next step is to keep solution within bounds
        x1 = (x1 .< l).*l+(l .<= x1).*(x1 .<= u).*x1+(u .< x1).*u			
        fx1 = problem.f(x1)

        status.f_calls += 1

        df  = fx1 - parameters.fx

        # If the function variation,df, is <0 we take test point as current
        # point. And if df>0 we use Metropolis condition to accept or
        # reject the test point as current point.
        if (df < 0 || rand() < exp(-T*df/(abs(parameters.fx)) / TolFun))
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

