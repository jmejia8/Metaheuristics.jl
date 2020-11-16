#_________________________________________________________________________#
#  Grasshopper Optimization Algorithm (GOA) source codes demo V0.1        #
#                                                                         #
#  Developed in Julia 0.6                                                 #
#                                                                         #
#  Programmer: Jesús Mejía                                                #
#                                                                         #
#  base on matlab code of S. Mirjalili                                    #
#                                                                         #
#  Main paper: S. Saremi, S. Mirjalili, A. Lewis                          #
#              Grasshopper Optimisation Algorithm: Theory and Application #
#               Advances in Engineering Software , in press,              #
#               DOI: http://dx.doi.org/10.1016/j.advengsoft.2017.01.004   #
#                                                                         #
#_________________________________________________________________________#

function initializationGOA(N,dim,up,down)

    if size(up,1)==1
        return rand(N, dim) .* (up - down) + down
    end

 
    X = zeros(N, dim)
    if size(up,1) > 1
        for i = 1:dim
            high    = up[i]
            low     = down[i]
            X[:, i] = rand(1,N) .* (high - low) + low
        end
    end

    return X
end

function S_func(r)
    f=0.5
    l=1.5
    return f* exp(-r / l) - exp(-r)  # Eq. (2.3) in the paper
end


# The Grasshopper Optimization Algorithm
function GOA(N, Max_evals, lb,ub, dim, fobj)

    Max_iter = div(Max_evals, N)
    # tic
    # println('GOA is now estimating the global optimum for your problem....')

    flag=0
    if size(ub,1)==1
        ub=ones(dim,1) * ub
        lb=ones(dim,1) * lb
    end

    if ( dim % 2  != 0) # this algorithm should be run with a even number of variables. This line is to handle odd number of variables
        dim = dim + 1
        ub = [ub 100]
        lb = [lb -100]
        flag=1
    end

    #Initialize the population of grasshoppers
    GrassHopperPositions = initializationGOA(N,dim,ub,lb)
    GrassHopperFitness   = zeros(N)

    fitness_history  = zeros(N,Max_iter)
    position_history = zeros(N,Max_iter,dim)
    Convergence_curve= zeros(1,Max_iter)
    Trajectories     = zeros(N,Max_iter)

    cMax = 1
    cMin = 0.00004
    #Calculate the fitness of initial grasshoppers

    for i = 1:size(GrassHopperPositions,1)
        if flag == 1
            GrassHopperFitness[i] = fobj(GrassHopperPositions[i,1:end-1])
        else
            GrassHopperFitness[i] = fobj(GrassHopperPositions[i, :])
        end

        fitness_history[i, 1]   = GrassHopperFitness[i]
        position_history[i, 1,:]= GrassHopperPositions[i, :]
        Trajectories[:, 1]      = GrassHopperPositions[:, 1]

    end

    sorted_indexes = sortperm(GrassHopperFitness)
    sorted_fitness = GrassHopperFitness[sorted_indexes]

    # Find the best grasshopper (target) in the first population 
    Sorted_grasshopper = zeros(N, dim)
    for newindex=1:N
        Sorted_grasshopper[newindex, :] = GrassHopperPositions[sorted_indexes[newindex],:]
    end

    TargetPosition = Sorted_grasshopper[1, :]
    TargetFitness  = sorted_fitness[1]

    # Main loop
    l = 2 # Start from the second iteration since the first iteration was dedicated to calculating the fitness of antlions
    while l < Max_iter + 1
        
        c = cMax-l*((cMax-cMin)/Max_iter) # Eq. (2.8) in the paper
        
         
         #################################################################################
         GrassHopperPositions_temp = zeros(N, dim)
          for i = 1:size(GrassHopperPositions,1)
            temp= GrassHopperPositions'
           # for k=1:2:dim  
                S_i = zeros(dim,1)
                for j = 1:N
                    if i != j
                        a_ = temp[:, j]
                        b_ = temp[:, i]

                        #println(size(ab_))

                        Dist = norm(a_ - b_) # Calculate the distance between two grasshoppers
                        
                        r_ij_vec = (a_ - b_) / (Dist)#+eps) # xj-xi/dij in Eq. (2.7)
                        xj_xi    = 2 + Dist % 2 # |xjd - xid| in Eq. (2.7) 
                        
                        s_ij = ((ub - lb)*c/2)*S_func(xj_xi) .* r_ij_vec # The first part inside the big bracket in Eq. (2.7)
                        S_i  = S_i+s_ij
                    end
                end
                S_i_total = S_i
                
          #  end
            
            X_new = c * S_i_total + (TargetPosition) # Eq. (2.7) in the paper      
            GrassHopperPositions_temp[i, :] = X_new' 
          end
          
        ########################################################################################
        # GrassHopperPositions
        GrassHopperPositions = GrassHopperPositions_temp
        
        for i=1:size(GrassHopperPositions,1)
            # Relocate grasshoppers that go outside the search space 
            # Tp = GrassHopperPositions[i, :] > ub'
            # Tm = GrassHopperPositions[i, :] < lb'
            # GrassHopperPositions[i, :] = (GrassHopperPositions[i, :] .* ( !( Tp+Tm)))+ub'.*Tp+lb'.*Tm
            
            # Calculating the objective values for all grasshoppers
            if flag == 1
                GrassHopperFitness[i] = fobj(GrassHopperPositions[i,1:end-1])
            else
                GrassHopperFitness[i] = fobj(GrassHopperPositions[i, :])
            end
            fitness_history[i, l]   = GrassHopperFitness[i]
            position_history[i, l,:]= GrassHopperPositions[i, :]
            
            Trajectories[:, l] = GrassHopperPositions[:, 1]
            
            # Update the target
            if GrassHopperFitness[i] < TargetFitness
                TargetPosition = GrassHopperPositions[i, :]
                TargetFitness  = GrassHopperFitness[i]
            end
        end
            
        Convergence_curve[l] = TargetFitness
        
        l = l + 1
    end


    if (flag==1)
        TargetPosition = TargetPosition[1:dim - 1]
    end

    # TargetFitness,TargetPosition,Convergence_curve,Trajectories,fitness_history, position_history
    return TargetFitness, TargetPosition, Convergence_curve
end