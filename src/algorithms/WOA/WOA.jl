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

mutable struct WOA
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
+=========== RESULT ==========+
| Iter.: 499
| f(x) = 4.56174e-104
| solution.x = [-1.04059445339676e-52, 1.7743142412652892e-52, 5.750781222647098e-53]
| f calls: 15000
| Total time: 0.0844 s
+============================+

julia> optimize(f, [-1 -1 -1; 1 1 1.0], WOA(N = 100))
+=========== RESULT ==========+
| Iter.: 499
| f(x) = 1.29795e-147
| solution.x = [1.306372696781744e-74, -3.017649118559932e-75, 3.3439182063846375e-74]
| f calls: 50000
| Total time: 0.1894 s
+============================+

```
"""
function WOA(;N = 30, information = Information(), options = Options())
  parameters = WOA(N)


  Algorithm(
    parameters,
    initialize! = initialize_woa!,
    update_state! = update_state_woa!,
    is_better = is_better,
    stop_criteria = stop_check,
    final_stage! = final_stage_woa!,
    information = information,
    options = options,
  )

end

function initialize_woa!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
   )

    lb, ub= problem.bounds[1,:], problem.bounds[2,:]

    #Initialize the positions of search agents
	max_it = 500
	options.iterations = options.iterations == 0 ? max_it : options.iterations
    options.f_calls_limit = options.f_calls_limit == 0 ?
                            options.iterations * parameters.N : options.f_calls_limit


    P = generate_population(problem.f, parameters.N, problem.bounds)
	status.population = P
    status.f_calls = parameters.N
    status.best_sol = deepcopy(get_best(status.population))
	
end

function update_state_woa!(
    problem,
    engine,
    parameters,
    status,
    information,
    options,
    iteration,
   )


  # Return back the search agents that go beyond the boundaries of the search space
  Leader_score = minimum(status)
  Leader_pos = minimizer(status)

  Max_iter = options.iterations
  N = parameters.N
  D = size(problem.bounds, 2)
  t = iteration

  a=2-t*((2)/Max_iter) # a decreases linearly fron 2 to 0 in Eq. (2.3)

  # a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
  a2=-1+t*((-1)/Max_iter)

  # Update the Position of search agents 
  for i=1:N
    r1 = rand() # r1 is a random number in [0,1]
    r2 = rand() # r2 is a random number in [0,1]

    A = 2a*r1 - a  # Eq. (2.3) in the paper
    C = 2r2        # Eq. (2.4) in the paper


    b = 1               #  parameters in Eq. (2.5)
    l = (a2-1)*rand()+1   #  parameters in Eq. (2.5)

    p = rand()        # p in Eq. (2.6)

    x = copy(status.population[i].x)
    for j = 1:D
      if p < 0.5   
        if abs(A) >= 1
          rand_leader_index = floor(Int, N* rand()+1)
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

    x = reset_to_violated_bounds!(x, problem.bounds)
    status.population[i] = generateChild(x, problem.f(x))
    status.f_calls += 1

    if engine.is_better(status.population[i], status.best_sol)
      status.best_sol = deepcopy(status.population[i])
    end
  end # for i

  status.stop = engine.stop_criteria(status, information, options) 

end

function final_stage_woa!(status, information, options)
    status.final_time = time()
end


# The Whale Optimization Algorithm

