"""
      optimize(
            f::Function, # objective function
            bounds::AbstractMatrix,
            method::AbstractAlgorithm = ECA();
            logger::Function = (status) -> nothing,
      )

Minimize a n-dimensional function `f` with domain `bounds` (2×n matrix) using `method = ECA()` by default.

# Example
Minimize f(x) = Σx² where x ∈ [-10, 10]³.

Solution:

```jldoctest
julia> f(x) = sum(x.^2)
f (generic function with 1 method)

julia> bounds = [  -10.0 -10 -10; # lower bounds
                    10.0  10 10 ] # upper bounds
2×3 Array{Float64,2}:
 -10.0  -10.0  -10.0
  10.0   10.0   10.0

julia> result = optimize(f, bounds)
+=========== RESULT ==========+
| Iter.: 1008
| f(x) = 6.48646e-163
| solution.x = [-4.054471688602619e-82, 4.2565448859996416e-82, 5.505242086898758e-82]
| f calls: 21187
| Total time: 0.1231 s
+============================+
```
"""

function optimize(
      f::Function, # objective function
      bounds::AbstractMatrix,
      method::AbstractAlgorithm = ECA();
      logger::Function = (status) -> nothing,
)

      problem = Problem(f, Array(bounds))
      convergence = State[]
      seed!(method.options.seed)
      method.options.debug && @info("Initializing population...")


      method.status.start_time = time()
      initialize!(
            method.status,
            method.parameters,
            problem,
            method.information,
            method.options
      )

      #####################################
      # common methods
      #####################################
      status = method.status
      information = method.information
      options = method.options
      ###################################

      ###################################
      # store convergence
      ###################################
      if options.store_convergence
            update_convergence!(convergence, status)
      end

      method.options.debug && @info("Starting main loop...")

      status.iteration = 0
      logger(status)

      while !status.stop
            status.iteration += 1

            update_state!(
                          method.status,
                          method.parameters,
                          problem,
                          method.information,
                          method.options
                         )

            if options.debug
                  status.final_time = time()
                  display(status)
            end

            if options.store_convergence
                  update_convergence!(convergence, status)
            end
            
            status.overall_time = time() - status.start_time
            logger(status)
            status.stop = status.stop || 
                        call_limit_stop_check(status, information, options) ||
                        iteration_stop_check(status, information, options)  ||
                        time_stop_check(status, information, options)
      end

      status.overall_time = time() - status.start_time

      final_stage!(
                   method.status,
                   method.parameters,
                   problem,
                   method.information,
                   method.options
                  )

      status.convergence = convergence

      return status

end

