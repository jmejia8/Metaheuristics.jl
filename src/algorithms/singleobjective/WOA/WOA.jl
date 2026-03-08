#_________________________________________________________________________#
#  Whale Optimization Algorithm (WOA) source codes demo 0.1               #
#                                                                         #
#                                                                         #
#  Programmer: Jesús Mejía                                                #
#  base on matlab code of S. Mirjalili                                    #
#                                                                         #
#                                                                         #
#   Main paper: S. Mirjalili, A. Lewis                                    #
#               The Whale Optimization Algorithm,                         #
#               Advances in Engineering Software , in press,              #
#               DOI: http://dx.doi.org/10.1016/j.advengsoft.2016.01.008   #
#                                                                         #
#_________________________________________________________________________#

mutable struct WOA <: AbstractParameters
  N::Int 
end

"""
    WOA(;N = 30, information = Information(), options = Options())

Parameters for the Whale Optimization Algorithm. `N` is the population size (number of
whales).


# Example


```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> optimize(f, [-1 -1 -1; 1 1 1.0], WOA())
Optimization Result
===================
  Iteration:       29
  Minimum:         3.05782e-62
  Minimizer:       [-2.47594e-32, -1.34764e-31, -1.08645e-31]
  Function calls:  870
  Total time:      0.0396 s
  Stop reason:     Due to Convergence Termination criterion.

julia> optimize(f, [-1 -1 -1; 1 1 1.0], WOA(N = 100))
Optimization Result
===================
  Iteration:       17
  Minimum:         1.38184e-66
  Minimizer:       [-1.12615e-33, -2.29611e-34, -2.46808e-34]
  Function calls:  1700
  Total time:      0.0015 s
  Stop reason:     Due to Convergence Termination criterion.

```
"""
function WOA(;N = 30, kargs...)
  parameters = WOA(N)


  Algorithm( parameters; kargs...)

end

function initialize!(
    status,
    parameters::WOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
  )

  lb = problem.search_space.lb
  ub = problem.search_space.ub

  #Initialize the positions of search agents
  max_it = 500
  options.iterations = options.iterations == 0 ? max_it : options.iterations
  options.f_calls_limit = options.f_calls_limit == 0 ?
  options.iterations * parameters.N : options.f_calls_limit


  # P = generate_population(parameters.N, problem)
  # best_sol = deepcopy(get_best(P))
  # status = State(best_sol, P)
  #=
  status.population = P
  status.best_sol = deepcopy(get_best(status.population))
  =#
  return gen_initial_state(problem,parameters,information,options,status)

end

function update_state!(
    status::State,
    parameters::WOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
  )


  # Return back the search agents that go beyond the boundaries of the search space
  Leader_score = minimum(status)
  Leader_pos = minimizer(status)

  Max_iter = options.iterations
  N = parameters.N
  D = getdim(problem)
  t = status.iteration
  rng = options.rng

  a=2-t*((2)/Max_iter) # a decreases linearly fron 2 to 0 in Eq. (2.3)

  # a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
  a2=-1+t*((-1)/Max_iter)

  X_new = zeros(N,D)

  # Update the Position of search agents 
  for i=1:N
    r1 = rand(rng) # r1 is a random number in [0,1]
    r2 = rand(rng) # r2 is a random number in [0,1]

    A = 2a*r1 - a  # Eq. (2.3) in the paper
    C = 2r2        # Eq. (2.4) in the paper


    b = 1               #  parameters in Eq. (2.5)
    l = (a2-1)*rand(rng)+1   #  parameters in Eq. (2.5)

    p = rand(rng)        # p in Eq. (2.6)

    x = zeros(D)
    for j = 1:D
      if p < 0.5   
        if abs(A) >= 1
          rand_leader_index = floor(Int, N* rand(rng)+1)
          X_rand  = status.population[rand_leader_index].x
          D_X_rand= abs(C*X_rand[j] - x[j]  ) # Eq. (2.7)
          x[j] = X_rand[j]-A*D_X_rand       # Eq. (2.8)

        else
          D_Leader=abs(C*Leader_pos[j] - x[j] ) # Eq. (2.1)
          x[j] =Leader_pos[j]-A*D_Leader      # Eq. (2.2)
        end

      else

        distance2Leader=abs(Leader_pos[j]-x[j])
        # Eq. (2.5)
        x[j] =distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos[j]

      end

    end # for j

    reset_to_violated_bounds!(x, problem.search_space)
    X_new[i,:] = x
  end # for i


  status.population = create_solutions(X_new, problem)
  sol = get_best(status.population)
  if is_better(sol, status.best_sol)
    status.best_sol = deepcopy(sol)
  end

end

function final_stage!(
    status::State,
    parameters::WOA,
    problem::AbstractProblem,
    information::Information,
    options::Options,
    args...;
    kargs...
  )
  status.final_time = time()
end


# The Whale Optimization Algorithm

