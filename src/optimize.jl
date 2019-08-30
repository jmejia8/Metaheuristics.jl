function optimize(f::Function, # objective function
                  bounds::Array,
                  method::AbstractAlgorithm
                  )

      problem = Problem(f, bounds)
      engine = method.engine

      method.options.debug && @info("Initializing population...")
      engine.initialize!(problem, engine, method.parameters, method.status, method.information, method.options)

      #####################################
      # common methods
      #####################################
      status = method.status
      information = method.information
      options     = method.options
      update_state! = engine.update_state!
      final_stage!  = engine.final_stage!
      ###################################

      ###################################
      # store convergence
      ###################################
      if options.store_convergence
            st = deepcopy(status)
            empty!(st.convergence)
            push!(status.convergence, st)
      end
      
      method.options.debug && @info("Starting main loop...")

      status.iteration = 0
      while !engine.stop_criteria(status, information, options)
            status.iteration += 1

            update_state!(problem,engine,method.parameters,method.status,method.information,method.options,status.iteration)
            
            options.debug && display(status)


            if options.store_convergence
                  st = deepcopy(status)
                  empty!(st.convergence)
                  push!(status.convergence, st)
            end
            
      end

      final_stage!(status, information, options)

      return status

end
