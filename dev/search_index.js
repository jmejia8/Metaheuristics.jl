var documenterSearchIndex = {"docs":
[{"location":"api/#API-References","page":"API References","title":"API References","text":"","category":"section"},{"location":"api/","page":"API References","title":"API References","text":"optimize","category":"page"},{"location":"api/#Metaheuristics.optimize","page":"API References","title":"Metaheuristics.optimize","text":"optimize(f::Function, bounds::Matrix{Float64}, method)\n\nMinimize a n-dimensional function f with domain bounds (2×n matrix) using method = ECA() by default.\n\nExample\n\nMinimize f(x) = Σx² where x ∈ [-10, 10]³.\n\nSolution:\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> bounds = [  -10.0 -10 -10; # lower bounds\n                    10.0  10 10 ] # upper bounds\n2×3 Array{Float64,2}:\n -10.0  -10.0  -10.0\n  10.0   10.0   10.0\n\njulia> result = optimize(f, bounds)\n+=========== RESULT ==========+\n| Iter.: 1008\n| f(x) = 6.48646e-163\n| solution.x = [-4.054471688602619e-82, 4.2565448859996416e-82, 5.505242086898758e-82]\n| f calls: 21187\n| Total time: 0.1231 s\n+============================+\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API References","title":"API References","text":"State","category":"page"},{"location":"api/#Metaheuristics.State","page":"API References","title":"Metaheuristics.State","text":"State datatype\n\nState is used to store the current metaheuristic status. In fact, the optimize function returns a State.\n\nbest_sol Stores the best solution found so far.\npopulation is an Array{typeof(best_sol)} for population-based algorithms.\nf_calls is the number of objective functions evaluations.\ng_calls  is the number of inequality constraints evaluations.\nh_calls is the number of equality constraints evaluations.\niteration is the current iteration.\nsuccess_rate percentage of new generated solutions better that their parents. \nconvergence used save the State at each iteration.\nstart_time saves the time() before the optimization proccess.\nfinal_time saves the time() after the optimization proccess.\nstop if true, then stops the optimization proccess.\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> bounds = [  -10.0 -10 -10; # lower bounds\n                    10.0  10 10 ] # upper bounds\n2×3 Array{Float64,2}:\n -10.0  -10.0  -10.0\n  10.0   10.0   10.0\n\njulia> state = optimize(f, bounds)\n+=========== RESULT ==========+\n| Iter.: 1009\n| f(x) = 7.16271e-163\n| solution.x = [-7.691251412064516e-83, 1.0826961235605951e-82, -8.358428300092186e-82]\n| f calls: 21190\n| Total time: 0.2526 s\n+============================+\n\njulia> minimum(state)\n7.162710802659093e-163\n\njulia> minimizer(state)\n3-element Array{Float64,1}:\n -7.691251412064516e-83\n  1.0826961235605951e-82\n -8.358428300092186e-82\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API References","title":"API References","text":"Information","category":"page"},{"location":"api/#Metaheuristics.Information","page":"API References","title":"Metaheuristics.Information","text":"Information Structure\n\nInformation can be used to store the true optimum in order to stop a metaheuristic early.\n\nProperties:\n\nf_optimum known minimum.\nx_optimum known minimizer.\n\nIf Options is provided, then optimize will stop when |f(x) - f(x_optimum)| < Options.f_tol or ‖ x - x_optimum ‖ < Options.x_tol (euclidean distance).\n\nExample\n\nIf you want an approximation to the minimum with accuracy of 1e-3 (|f(x) - f(x*)| < 1e-3), then you may use Information.\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> bounds = [  -10.0 -10 -10; # lower bounds\n                    10.0  10 10 ] # upper bounds\n2×3 Array{Float64,2}:\n -10.0  -10.0  -10.0\n  10.0   10.0   10.0\n\njulia> information = Information(f_optimum = 0.0)\nInformation(0.0, Float64[])\n\njulia> options = Options(f_tol = 1e-3)\nOptions(0.0, 0.001, 0.0, 0.0, 1000.0, 0.0, 0.0, 0, false, true, false, :minimize)\n\njulia> state = optimize(f, bounds, ECA(information=information, options=options))\n+=========== RESULT ==========+\n| Iter.: 22\n| f(x) = 0.000650243\n| solution.x = [0.022811671589729583, 0.007052331140376011, -0.008951836265056107]\n| f calls: 474\n| Total time: 0.0106 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API References","title":"API References","text":"Options","category":"page"},{"location":"api/#Metaheuristics.Options","page":"API References","title":"Metaheuristics.Options","text":"Options()\n\nOptions stores common settings for metaheuristics such as the maximum number of iterations debug options, maximum number of function evaluations, etc.\n\nMain properties:\n\nx_tol tolerance to the true minimizer if specified in Information.\nf_tol tolerance to the true minimum if specified in Information.\nf_calls_limit is the maximum number of function evaluations limit.\niterations is the maximum number iterationn permited.\nstore_convergence if true, then push the current State in State.convergence at each generation/iteration\ndebug if true, then optimize function reports the current State (and interest information) for each iterations.\n\nExample\n\njulia> options = Options(f_calls_limit = 1000, debug=true)\nOptions(0.0, 0.0, 0.0, 0.0, 1000.0, 0.0, 0.0, 0, false, true, true, :minimize)\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> bounds = [  -10.0 -10 -10; # lower bounds\n                    10.0  10 10 ] # upper bounds\n2×3 Array{Float64,2}:\n -10.0  -10.0  -10.0\n  10.0   10.0   10.0\n\njulia> state = optimize(f, bounds, ECA(options=options))\n[ Info: Initializing population...\n[ Info: Starting main loop...\n+=========== RESULT ==========+\n| Iter.: 1\n| f(x) = 17.4821\n| solution.x = [-0.21198067919245656, -4.0408459028606885, -1.0529741414392084]\n| f calls: 42\n| Total time: 0.0003 s\n+============================+\n\n...\n\n[ Info: Stopped since call_limit was met.\n+=========== RESULT ==========+\n| Iter.: 47\n| f(x) = 1.0924e-06\n| solution.x = [-0.0008050045939313401, 0.000319255968803667, -0.0005851867103384286]\n| f calls: 1000\n| Total time: 0.0258 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"api/","page":"API References","title":"API References","text":"convergence","category":"page"},{"location":"api/#Metaheuristics.convergence","page":"API References","title":"Metaheuristics.convergence","text":"convergence(state)\n\nget the data (touple with the number of function evaluations and fuction values) to plot the convergence graph. \n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> bounds = [  -10.0 -10 -10; # lower bounds\n                    10.0  10 10 ] # upper bounds\n2×3 Array{Float64,2}:\n -10.0  -10.0  -10.0\n  10.0   10.0   10.0\n\njulia> state = optimize(f, bounds, ECA(options=Options(store_convergence=true)))\n+=========== RESULT ==========+\n| Iter.: 1022\n| f(x) = 7.95324e-163\n| solution.x = [-7.782044850211721e-82, 3.590044165897827e-82, -2.4665318114710003e-82]\n| f calls: 21469\n| Total time: 0.3300 s\n+============================+\n\njulia> n_fes, fxs = convergence(state);\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API References","title":"API References","text":"minimizer","category":"page"},{"location":"api/#Metaheuristics.minimizer","page":"API References","title":"Metaheuristics.minimizer","text":"minimizer(state)\n\nReturns the approximation to the minimizer (argmin f(x)) stored in state.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API References","title":"API References","text":"minimum(state::State)","category":"page"},{"location":"api/#Base.minimum-Tuple{State}","page":"API References","title":"Base.minimum","text":"minimum(state::Metaheuristics.State)\n\nReturns the approximation to the minimum (min f(x)) stored in state.\n\n\n\n\n\n","category":"method"},{"location":"api/","page":"API References","title":"API References","text":"positions","category":"page"},{"location":"api/#Metaheuristics.positions","page":"API References","title":"Metaheuristics.positions","text":"positions(state)\n\nIf state.population has N solutions, then returns a N×d Matrix.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API References","title":"API References","text":"fvals","category":"page"},{"location":"api/#Metaheuristics.fvals","page":"API References","title":"Metaheuristics.fvals","text":"fvals(state)\n\nIf state.population has N solutions, then returns a Vector with the  objective function values from items in state.population.\n\n\n\n\n\n","category":"function"},{"location":"api/","page":"API References","title":"API References","text":"nfes","category":"page"},{"location":"api/#Metaheuristics.nfes","page":"API References","title":"Metaheuristics.nfes","text":"nfes(state)\n\nget the number of function evaluations.\n\n\n\n\n\n","category":"function"},{"location":"#Introduction","page":"Introduction","title":"Introduction","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Optimization is one of the most common task in the scientific and industry field but real-world problems require high-performance algorithms to optimize non-differentiable, non-convex, dicontinuous functions. Different metaheuristics algorithms have been proposed to solve optimization problems but without strong assumptions about the objective function.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This package implements state-of-the-art metaheuristics algorithms for global optimization. The aim of this package is to provide easy to use (and fast) metaheuristics for numerical global optimization.","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Open the Julia (Julia 0.7 or Later) REPL and press ] to open the Pkg prompt. To add this package, use the add command:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"pkg> add https://github.com/jmejia8/Metaheuristics.jl.git","category":"page"},{"location":"#Quick-Start","page":"Introduction","title":"Quick Start","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Assume you want to solve the following minimization problem.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: Rastrigin Surface)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Minimize:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"f(x) = 10D + sum_i=1^D  x_i^2 - 10cos(2pi x_i)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"where xin-5 5^D, i.e., -5 leq x_i leq 5 for i=1ldotsD. D is the dimension number, assume D=10.","category":"page"},{"location":"#Solution","page":"Introduction","title":"Solution","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Firstly, import the Metaheuristics package:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"using Metaheuristics","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Code the objective function:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x)  )","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Instantiate the bounds, note that bounds should be a 2times 10 Matrix where the first row corresponds to the lower bounds whilst the second row corresponds to the upper bounds.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"D = 10\nbounds = [-5ones(D) 5ones(D)]'","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Approximate the optimum using the function optimize.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"result = optimize(f, bounds)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Optimize returns a State datatype which contains some information about the approximation. For instance, you may use mainly two functions to obtain such approximation.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"@show minimum(result)\n@show minimizer(result)","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"After reading this tutorial you'll become an expert using Metaheuristics module.","category":"page"},{"location":"tutorial/#Minimization-Problem","page":"Tutorial","title":"Minimization Problem","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Assume you want to optimize the following minimization problem:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Minimize:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"f(x) = 10D + sum_i=1^D x_i^2 - 10cos(2pi x_i)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where xin -5 5^D, that is, each coordinate in x is between -5 and 5. Use D=10.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Note that the global optimum is obtained when x_i = 0 for all i. Thus, min f(x) = 0.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Objective function:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"f(x) = 10length(x) + sum( x.^2 - 10cos.(2pi*x) )","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Bounds:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"bounds = [-5ones(10) 5ones(10)]'","category":"page"},{"location":"tutorial/#Providing-Information","page":"Tutorial","title":"Providing Information","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Since the optimum is known, then we can provide this information to the optimizer.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"information = Information(f_optimum = 0.0)","category":"page"},{"location":"tutorial/#Common-Settings","page":"Tutorial","title":"Common Settings","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Usually users could require to limit the number of generation/iteration or the number of function evaluations. To do that, let's assume that the metaheuristic should evaluate at most 9000D times the objective function. Moreover, since information is provided, then we can set the desired accuracy (f(x) - f(x^*) ) to 10^-5.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)","category":"page"},{"location":"tutorial/#Choose-a-Metaheuristic","page":"Tutorial","title":"Choose a Metaheuristic","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Metaheuristics.jl provides different metaheuristics for optimization such as Evolutionary Centers Algorithm (ECA), Differential Evolution (DE), Particle Swarm Optimization (PSO), etc. In this tutorial we will use ECA, but you can use another algorithm following but the same steps.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The metaheuristics accept its parameters but share two common and optional settings information and options.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"algorithm = ECA(information = information, options = options)","category":"page"},{"location":"tutorial/#Optimize","page":"Tutorial","title":"Optimize","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now, we are able to approximate the optimum. To do that is necessary to use the optimize function as follows:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"result = optimize(f, bounds, algorithm)","category":"page"},{"location":"tutorial/#Get-the-Results","page":"Tutorial","title":"Get the Results","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Once optimize stopped, then we can get the approximate solutions.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Approximated minimum:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"fx = minimum(result)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Approximated minimizer:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"x = minimizer(result)","category":"page"},{"location":"tutorial/#Get-Information-about-the-Resulting-Population","page":"Tutorial","title":"Get Information about the Resulting Population","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Sometimes is useful to analyze the resulting population (for population-based metaheuristics). To do that you can use fvals to get objective function evaluation and positions to get their positions.","category":"page"},{"location":"tutorial/#Bonus","page":"Tutorial","title":"Bonus","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"We recommend you to save  your program in a function for performance purposes:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using Metaheuristics\n\nfunction main()\n    # objective function\n    f(x) = 10length(x) + sum( x.^2 - 10cos.(2π*x) )\n    \n    # limits/bounds\n    bounds = [-5ones(10) 5ones(10)]'\n    \n    # information on the minimization problem\n    information = Information(f_optimum = 0.0)\n\n    # generic settings\n    options = Options(f_calls_limit = 9000*10, f_tol = 1e-5)\n    \n    # metaheuristic used to optimize\n    algorithm = ECA(information = information, options = options)\n\n    # start the minimization proccess\n    result = optimize(f, bounds, algorithm)\n\n    \n    fx = minimum(result)\n    x = minimizer(result)\n\n    @show fx\n    @show x\nend\n","category":"page"},{"location":"tutorial/#Summary","page":"Tutorial","title":"Summary","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now you are able to approximate global optimum solutions using Metaheuristics.","category":"page"},{"location":"algorithms/#Algorithms","page":"Algorithms","title":"Algorithms","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"List of implemented metaheuristics.","category":"page"},{"location":"algorithms/#Evolutionary-Centers-Algorithm","page":"Algorithms","title":"Evolutionary Centers Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"ECA","category":"page"},{"location":"algorithms/#Metaheuristics.ECA","page":"Algorithms","title":"Metaheuristics.ECA","text":"ECA(;\n    η_max = 2.0,\n    K = 7,\n    N = 0,\n    N_init = N,\n    p_exploit = 0.95,\n    p_bin = 0.02,\n    ε = 0.0,\n    p_cr = Float64[],\n    adaptive = false,\n    resize_population = false,\n    information = Information(),\n    options = Options()\n)\n\nParameters for the metaheuristic ECA: step-size η_max,K is number of vectors to generate the center of mass, N is the population size.\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], ECA())\n\n+=========== RESULT ==========+\n| Iter.: 1021\n| f(x) = 1.68681e-163\n| solution.x = [2.5517634463667404e-82, -2.9182760041942484e-82, -1.3565584801935802e-82]\n| f calls: 21454\n| Total time: 0.0894 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], ECA(N = 10, η_max = 1.0, K = 3))\n+=========== RESULT ==========+\n| Iter.: 1506\n| f(x) = 0.000172391\n| solution.x = [-6.340714627875324e-5, -0.004127226953894587, 0.012464071313908906]\n| f calls: 15069\n| Total time: 0.0531 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Differential-Evolution","page":"Algorithms","title":"Differential Evolution","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"DE","category":"page"},{"location":"algorithms/#Metaheuristics.DE","page":"Algorithms","title":"Metaheuristics.DE","text":"DE(;\n    N  = 0,\n    F  = 1.0,\n    CR = 0.9,\n    strategy = :rand1,\n    information = Information(),\n    options = Options()\n)\n\nParameters for Differential Evolution (DE) algorithm: step-size F,CR controlls the binomial crossover, N is the population size. The parameter trategy is related to the variation operator (:rand1, :rand2, :best1, :best2, :randToBest1).\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], DE())\n+=========== RESULT ==========+\n| Iter.: 437\n| f(x) = 0\n| solution.x = [0.0, 0.0, 0.0]\n| f calls: 13134\n| Total time: 0.3102 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], DE(N=50, F=1.5, CR=0.8))\n+=========== RESULT ==========+\n| Iter.: 599\n| f(x) = 9.02214e-25\n| solution.x = [-4.1003250484858545e-13, -6.090890160928905e-13, -6.025762626763004e-13]\n| f calls: 30000\n| Total time: 0.0616 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Particle-Swarm-Optimization","page":"Algorithms","title":"Particle Swarm Optimization","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"PSO","category":"page"},{"location":"algorithms/#Metaheuristics.PSO","page":"Algorithms","title":"Metaheuristics.PSO","text":"PSO(;\n    N  = 0,\n    C1 = 2.0,\n    C2 = 2.0,\n    ω  = 0.8,\n    information = Information(),\n    options = Options()\n)\n\nParameters for Particle Swarm Optimization (PSO) algorithm: learning rates C1 and C2, N is the population size and ω controlls the inertia weight. \n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO())\n+=========== RESULT ==========+\n| Iter.: 999\n| f(x) = 3.23944e-48\n| solution.x = [1.0698542573895642e-24, -1.4298101555926563e-24, -2.247029420442994e-25]\n| f calls: 30000\n| Total time: 0.4973 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], PSO(N = 100, C1=1.5, C2=1.5, ω = 0.7))\n+=========== RESULT ==========+\n| Iter.: 299\n| f(x) = 1.41505e-38\n| solution.x = [2.161357427851024e-20, -1.1599444038307776e-19, 1.5122345732802047e-20]\n| f calls: 30000\n| Total time: 0.2128 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Artificial-Bee-Colony","page":"Algorithms","title":"Artificial Bee Colony","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"ABC","category":"page"},{"location":"algorithms/#Metaheuristics.ABC","page":"Algorithms","title":"Metaheuristics.ABC","text":"ABC(;\n    N = 50,\n    Ne = div(N+1, 2),\n    No = div(N+1, 2),\n    limit=10,\n    information = Information(),\n    options = Options()\n)\n\nABC implements the original parameters for the Artificial Bee Colony Algorithm. N is the population size, Ne is the number of employees, No is the number of outlookers bees. limit is related to the times that a solution is visited.\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], ABC())\n+=========== RESULT ==========+\n| Iter.: 593\n| f(x) = 3.54833e-25\n| solution.x = [3.448700205761237e-13, 4.805851037329074e-13, 7.025504722610375e-14]\n| f calls: 30019\n| Total time: 0.2323 s\n+============================+\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], ABC(N = 80,  No = 20, Ne = 50, limit=5))\n+=========== RESULT ==========+\n| Iter.: 405\n| f(x) = 2.24846e-07\n| solution.x = [0.0002682351072804559, 0.00020460896416511776, 0.0003332131896109299]\n| f calls: 30043\n| Total time: 0.2652 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#MOEA/D-DE","page":"Algorithms","title":"MOEA/D-DE","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"MOEAD_DE","category":"page"},{"location":"algorithms/#Metaheuristics.MOEAD_DE","page":"Algorithms","title":"Metaheuristics.MOEAD_DE","text":"MOEAD_DE(D::Int, nobjectives::Int)\n\nMOEAD_DE implements the original version of MOEA/D-DE. It uses the contraint handling method based on the sum of violations (for constrained optimizaton): g(x, λ, z) = max(λ .* abs.(fx - z)) + sum(max.(0, gx)) + sum(abs.(hx))\n\nTo use MOEAD_DE, the output from the objective function should be a 3-touple (f::Vector, g::Vector, h::Vector), where f contains the objective functions, g and h are the equality and inequality constraints respectively.\n\nA feasible solution is such that g_i(x) ≤ 0 and h_j(x) = 0.\n\nExample\n\nAssume you want to solve the following optimizaton problem:\n\nMinimize:\n\nf(x) = (x_1, x_2)\n\nsubject to:\n\ng(x) = x_1^2 + x_2^2 - 1 ≤ 0\n\nx_1, x_2 ∈ [-1, 1]\n\nA solution can be:\n\n\n# Dimension\nD = 2\n\n# Objective function\nf(x) = ( x, [sum(x.^2) - 1], [0.0] ) \n\n# bounds\nbounds = [-1 -1;\n           1  1.0\n        ]\n\n# define the parameters\nmoead_de = MOEAD_DE(D, 2, N = 300, options=Options(debug=false, iterations = 500))\n\n# optimize\nstatus_moead = optimize(f, bounds, moead_de)\n\n# show results\ndisplay(status_moead)\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Gravitational-Search-Algorithm","page":"Algorithms","title":"Gravitational Search Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"CGSA","category":"page"},{"location":"algorithms/#Metaheuristics.CGSA","page":"Algorithms","title":"Metaheuristics.CGSA","text":"CGSA(;\n    N::Int    = 30,\n    chValueInitial::Real   = 20,\n    chaosIndex::Real   = 9,\n    ElitistCheck::Int    = 1,\n    Rpower::Int    = 1,\n    Rnorm::Int    = 2,\n    wMax::Real   = chValueInitial,\n    wMin::Real   = 1e-10,\n    information = Information(),\n    options = Options()\n)\n\nCGSA is an extension of the GSA algorithm but with Chaotic gravitational constants for the gravitational search algorithm.\n\nRef. Chaotic gravitational constants for the gravitational search algorithm. Applied Soft Computing 53 (2017): 407-419.\n\nParameters:\n\nN: Population size\nchValueInitial: Initial value for the chaos value\nchaosIndex: Integer 1 ≤ chaosIndex ≤ 10 is the function that model the chaos\nRpower: power related to the distance norm(x)^Rpower\nRnorm: is the value as in norm(x, Rnorm)\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], CGSA())\n+=========== RESULT ==========+\n| Iter.: 499\n| f(x) = 0.000235956\n| solution.x = [0.0028549782101697785, -0.0031385153631797724, 0.014763299731686608]\n| f calls: 15000\n| Total time: 0.1003 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], CGSA(N = 80, chaosIndex = 1))\n+=========== RESULT ==========+\n| Iter.: 499\n| f(x) = 0.000102054\n| solution.x = [0.00559987302269564, 0.00017535321765604905, 0.008406213942044265]\n| f calls: 40000\n| Total time: 0.5461 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Simulated-Annealing","page":"Algorithms","title":"Simulated Annealing","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"SA","category":"page"},{"location":"algorithms/#Metaheuristics.SA","page":"Algorithms","title":"Metaheuristics.SA","text":"    SA(;\n        x_initial::Vector = zeros(0),\n        N::Int = 500,\n        tol_fun::Real= 1e-4,\n        information = Information(),\n        options = Options()\n    )\n\nParameters for the method of Simulated Annealing (Kirkpatrick et al., 1983).\n\nParameters:\n\nx_intial: Inital solution. If empty, then SA will generate a random one within the bounds.\nN: The number of test points per iteration.\ntol_fun: tolerance value for the Metropolis condition to accept or reject the test point as current point.\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], SA())\n+=========== RESULT ==========+\n| Iter.: 60\n| f(x) = 2.84574e-73\n| solution.x = [-5.307880224731971e-37, -5.183298967486749e-38, 1.2301984439451926e-38]\n| f calls: 29502\n| Total time: 0.0465 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], SA(N = 100, x_initial = [1, 0.5, -1]))\n+=========== RESULT ==========+\n| Iter.: 300\n| f(x) = 1.29349e-70\n| solution.x = [-7.62307964668667e-36, 8.432089040013441e-36, -3.7077496015659554e-37]\n| f calls: 29902\n| Total time: 0.0466 s\n+============================+\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#Whale-Optimization-Algorithm","page":"Algorithms","title":"Whale Optimization Algorithm","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"WOA","category":"page"},{"location":"algorithms/#Metaheuristics.WOA","page":"Algorithms","title":"Metaheuristics.WOA","text":"WOA(;N = 30, information = Information(), options = Options())\n\nParameters for the Whale Optimization Algorithm. N is the population size (number of whales).\n\nExample\n\njulia> f(x) = sum(x.^2)\nf (generic function with 1 method)\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], WOA())\n+=========== RESULT ==========+\n| Iter.: 499\n| f(x) = 4.56174e-104\n| solution.x = [-1.04059445339676e-52, 1.7743142412652892e-52, 5.750781222647098e-53]\n| f calls: 15000\n| Total time: 0.0844 s\n+============================+\n\njulia> optimize(f, [-1 -1 -1; 1 1 1.0], WOA(N = 100))\n+=========== RESULT ==========+\n| Iter.: 499\n| f(x) = 1.29795e-147\n| solution.x = [1.306372696781744e-74, -3.017649118559932e-75, 3.3439182063846375e-74]\n| f calls: 50000\n| Total time: 0.1894 s\n+============================+\n\n\n\n\n\n\n","category":"type"},{"location":"algorithms/#NSGA-II","page":"Algorithms","title":"NSGA-II","text":"","category":"section"},{"location":"algorithms/","page":"Algorithms","title":"Algorithms","text":"NSGA2","category":"page"},{"location":"algorithms/#Metaheuristics.NSGA2","page":"Algorithms","title":"Metaheuristics.NSGA2","text":"function NSGA2(;\n    N = 100,\n    η_cr = 20,\n    p_cr = 0.9,\n    η_m = 20,\n    p_m = 1.0 / D,\n    ε = eps(),\n    information = Information(),\n    options = Options(),\n)\n\nParameters for the metaheuristic NSGA-II.\n\nParameters:\n\nN Population size.\nη_cr  η for the crossover.\np_cr Crossover probability.\nη_m  η for the mutation operator.\np_m Mutation probability (1/D for D-dimensional problem by default).\n\nTo use NSGA2, the output from the objective function should be a 3-touple (f::Vector, g::Vector, h::Vector), where f contains the objective functions, g and h are the equality and inequality constraints respectively.\n\nA feasible solution is such that g_i(x) ≤ 0 and h_j(x) = 0.\n\nusing Metaheuristics\n\n# Dimension\nD = 2\n\n# Objective function\nf(x) = ( x, [sum(x.^2) - 1], [0.0] ) \n\n# bounds\nbounds = [-1 -1;\n           1  1.0\n        ]\n\n# define the parameters (use `NSGA2()` for using default parameters)\nnsga2 = NSGA2(N = 100, p_cr = 0.85)\n\n# optimize\nstatus = optimize(f, bounds, nsga2)\n\n# show results\ndisplay(status)\n\n\n\n\n\n","category":"type"}]
}
