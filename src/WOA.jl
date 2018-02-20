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

# This function initialize the first population of search agents
function initialization(SearchAgents_no,dim,ub,lb)
    Boundary_no= size(ub,2) # numnber of boundaries

    # If the boundaries of all variables are equal and user enter a signle
    # number for both ub and lb
    if Boundary_no==1
        Positions=rand(SearchAgents_no,dim).*(ub-lb)+lb
    end

    # If each variable has a different lb and ub
    if Boundary_no > 1
        for i=1:dim
            ub_i = ub[i]
            lb_i = lb[i]
            Positions(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i
        end
    end
    return Positions 
end

# The Whale Optimization Algorithm
function WOA(SearchAgents_no, Max_evals, lb, ub, dim, fobj)

    Max_iter = div(Max_evals, SearchAgents_no)

    # initialize position vector and score for the leader
    Leader_pos  = zeros(1,dim)
    Leader_score= Inf #change this to -inf for maximization problems


    #Initialize the positions of search agents
    Positions = initialization(SearchAgents_no,dim,ub,lb)

    Convergence_curve = zeros(1,Max_iter)

    t=0# Loop counter

    # Main loop
    while t < Max_iter
        for i=1:size(Positions,1)
            
            # Return back the search agents that go beyond the boundaries of the search space
            # Flag4ub = Positions[i,:] .> ub
            # Flag4lb = Positions[i,:] .< lb
            # Positions[i,:]=( Positions[i,:] .* (!(Flag4ub + Flag4lb)) ) + ub .* Flag4ub + lb .* Flag4lb
            
            # Calculate objective function for each search agent
            fitness=fobj(Positions[i,:])
            
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
        for i=1:size(Positions,1)
            r1 = rand() # r1 is a random number in [0,1]
            r2 = rand() # r2 is a random number in [0,1]
            
            A = 2a*r1 - a  # Eq. (2.3) in the paper
            C = 2r2        # Eq. (2.4) in the paper
            
            
            b = 1               #  parameters in Eq. (2.5)
            l = (a2-1)*rand()+1   #  parameters in Eq. (2.5)
            
            p = rand()        # p in Eq. (2.6)
            
            for j = 1:size(Positions,2)
                
                if p < 0.5   
                    if abs(A) >= 1
                        rand_leader_index = floor(Integer, SearchAgents_no*rand()+1)
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
        Convergence_curve[t]=Leader_score

    end
    
    return Leader_score,Leader_pos,Convergence_curve


end


