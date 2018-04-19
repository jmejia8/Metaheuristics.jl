#_________________________________________________________________________#
#  Whale Optimization Algorithm (WOA) source codes demo 0.1               #
#                                                                         #
#  Developed in Julia(0.6)                                                #
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

# The Whale Optimization Algorithm
function WOA(fobj::Function,
                D::Int;
                N::Int = 30,
        max_evals::Int = 10000D, 
      showResults::Bool = true,
         saveLast::String = "",
  saveConvergence::String = "",
             limits = [-100.0, 100.0])

    Max_iter = div(max_evals, N) + 1
    lb, ub = limits

    # bounds vectors
    lb, ub = limits[1,:], limits[2,:]
    if length(lb) < D
        lb = ones(D) * lb[1]
        ub = ones(D) * ub[1]
    end


    # initialize position vector and score for the leader
    Leader_pos  = zeros(1,D)
    Leader_score= Inf #change this to -inf for maximization problems


    #Initialize the positions of search agents
    Positions = initializePop(N, D, lb, ub)

    convergence = []

    t=0# Loop counter
    nevals = 0

    # Main loop
    while t < Max_iter
        # Return back the search agents that go beyond the boundaries of the search space
        Positions = correctPop(Positions, lb, ub)

        for i=1:N
            

            # Calculate objective function for each search agent
            fitness=fobj(Positions[i,:])
            nevals += 1
            
            # Update the leader
            if fitness < Leader_score # Change this to > for maximization problem
                Leader_score = fitness # Update alpha
                Leader_pos   = Positions[i,:]
            end
            
        end
        
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
            
            for j = 1:D
                
                if p < 0.5   
                    if abs(A) >= 1
                        rand_leader_index = floor(Integer, N* rand()+1)
                        X_rand  = Positions[rand_leader_index, :]
                        D_X_rand= abs(C*X_rand[j]-Positions[i,j]) # Eq. (2.7)
                        Positions[i,j]=X_rand[j]-A*D_X_rand       # Eq. (2.8)
                        
                    elseif abs(A)<1
                        D_Leader=abs(C*Leader_pos[j]-Positions[i,j]) # Eq. (2.1)
                        Positions[i,j]=Leader_pos[j]-A*D_Leader      # Eq. (2.2)
                    end
                    
                elseif p>=0.5
                  
                    distance2Leader=abs(Leader_pos[j]-Positions[i,j])
                    # Eq. (2.5)
                    Positions[i,j]=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos[j]
                    
                end
                
            end
        end
        t=t+1

        push!(convergence, [nevals Leader_score])

    end

    if saveConvergence != ""
       writecsv(saveConvergence, convergence)
    end

    if saveLast != ""
       writecsv(saveLast, Positions)
    end
    
    if showResults
        println("===========[ ECA results ]=============")
        println("| Generations = $t")
        println("| Evals       = ", t*N)
        @printf("| best f.     = %e\n", Leader_score)
        println("=======================================")
    end

    return Leader_pos, Leader_score


end


