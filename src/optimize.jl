function optimize(f::Function, # objective function
                  bounds::Array,
                  method::AbstractAlgorithm
                  )

      problem = Problem(f, bounds)
      engine = method.engine
      convergence = State[]
      method.options.debug && @info("Initializing population...")

      method.status.start_time = time()
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
            update_convergence!(convergence, status)
      end
      
      method.options.debug && @info("Starting main loop...")

      status.iteration = 0
      while !status.stop
            status.iteration += 1

            update_state!(problem,engine,method.parameters,method.status,method.information,method.options,status.iteration)
            
            if options.debug
                  status.final_time = time()
                  display(status)
            end

            if options.store_convergence
                  update_convergence!(convergence, status)
            end
            
      end

      final_stage!(status, information, options)

      status.convergence = convergence

      return status

end
