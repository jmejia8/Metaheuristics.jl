function Cal_par(prob_k)
    # CEC2021 Real world multi-objective Constrained Optimization Test Suite 
    # Abhishek Kumar (email: abhishek.kumar.eee13@iitbhu.ac.in, Indian Institute of Technology (BHU), Varanasi) 

    # prob_k -> Index of problem.
    # n  -> Dimension of the problem.
    # fn -> Number of objective.
    # g  -> Number of inequility constraints.
    # h  -> Number of equality constraints.
    # xmin -> lower bound of decision variables.
    # xmax -> upper bound of decision variables.


    D        = [4,5,3,4,4,7,4,7,4,2,3,4,7,5,3,2,6,3,10,4,6,9,6,9,2,3,3,7,7,25,25,25,30,30,30,28,28,28,28,34,34,34,34,34,34,34,18,18,18,6];
    O        = [2,2,2,2,2,2,2,3,2,2,5,2,3,2,2,2,3,2,3,2,2,2,2,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,2,3,2,2,3,3,4,2,2,3,2];
    fn   = O[prob_k];
    n    = D[prob_k]; 
    gn       = [2,5,3,4,4,11,1,9,1,2,7,1,11,8,8,2,9,3,10,7,4,2,1,0,2,1,3,4,9,24,24,24,29,29,29,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
    hn       = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,4,6,0,1,0,4,0,0,0,0,0,0,0,24,24,24,24,26,26,26,26,26,26,26,12,12,12,1];
    g     = gn[prob_k];
    h     = hn[prob_k];
    ## range
    # bound constraint definitions for all 18 test functions
    xmin  = [ zeros(1) for i in 1:50 ]
    xmax  = [ zeros(1) for i in 1:50 ]
    xmin[1]    = [0.51,0.51,10,10];
    xmax[1]    = [99.49,99.49,200,200];
    xmin[2]    = [0.05,0.2,0.2,0.35,3];
    xmax[2]    = [0.5,0.5,0.6,0.5,6];
    xmin[3]    = [1e-5,1e-5,1];
    xmax[3]    = [100, 100, 3];
    xmin[4]    = [0.125,0.1,0.1,0.125];
    xmax[4]    = [5,10,10,5];
    xmin[5]    = [55,75,1000,11];
    xmax[5]    = [80,110,3000,20];
    xmin[6]    = [2.6,0.7,16.51,7.3,7.3,2.9,5];
    xmax[6]    = [3.6,0.8,28.49,8.3,8.3,3.9,5.5];
    xmin[7]    = [11.51,11.51,11.51,11.51];
    xmax[7]    = [60.49,60.49,60.49,60.49];
    xmin[8]    = [0.5,0.45,0.5,0.5,0.875,0.4,0.4];
    xmax[8]    = [1.5,1.35,1.5,1.5,2.625,1.2,1.2];
    xmin[9]    = [1,sqrt(2),sqrt(2),1];
    xmax[9]    = [3,3,3,3];
    xmin[10]   = [0.1,0.5];
    xmax[10]   = [2,2.5];
    xmin[11]   = [0.01,0.01, 0.01];
    xmax[11]   = [0.45,0.1, 0.1];
    xmin[12]   = [10, 10, 0.9, 0.9];
    xmax[12]   = [80, 50, 5, 5];  
    xmin[13]    = [2.6,0.7,16.51,7.3,7.3,2.9,5];
    xmax[13]    = [3.6,0.8,28.49,8.3,8.3,3.9,5.5];
    xmin[14]    =  [60,90, 1, 0, 2];
    xmax[14]    = [80, 110, 3, 1000, 9];
    xmin[15]    = [0.51,0.6,0.51];
    xmax[15]    = [70.49,3,42.49];
    xmin[16]    = [0.01, 0.20];
    xmax[16]    = [0.05, 1];
    xmin[17]    = [150.0, 20.0, 13.0, 10.0, 14.0, 0.63];
    xmax[17]    = [274.32, 32.31, 25.0, 11.71, 18.0, 0.75];
    xmin[18]    = [136,56,1.4];
    xmax[18]    = [146,68,2.2];
    xmin[19]   = [ 0.51,0.51,0.51,250,250,250,6,4,40,10];
    xmax[19]   = [3.49,3.49,3.49,2500,2500,2500,20,16,700,450];
    xmin[20]   = [ 1, 1,  1e-6,1];
    xmax[20]   = [16, 16, 16*1e-6,16];
    xmin[21] = [1.3, 2.5, 1.3, 1.3, 1.3, 1.3];
    xmax[21] = [1.7, 3.5, 1.7, 1.7, 1.7, 1.7];
    xmin[22] = zeros(9);
    xmax[22] = [100,200,100,100,100,100,200,100,200];
    xmin[23] = [0,0,0,0,0.00001,0.00001];
    xmax[23] = [ 1,1,1,1,16,16];
    xmin[24] = [0,0,0,0,1000,0,100,100,100];
    xmax[24] = [10,200,100,200,2000000,600,600,600,900];
    xmin[25] = [0,-0.49];
    xmax[25] = [1.6,1.49];
    xmin[26] = [0.5,0.5,-0.49];
    xmax[26] = [1.4,1.4,1.49];
    xmin[27] = [0.2,-2.22554,-0.49];
    xmax[27] = [1,-1,1.49];
    xmin[28] = [0,0,0,0,-0.49,-0.49,0];
    xmax[28] = [20,20,10,10,1.49,1.49,40];
    xmin[29] = [0,0,0,-0.49,-0.49,-0.49,-0.49];
    xmax[29] = [100,100,100,1.49,1.49,1.49,1.49];
    xmin[30]   = -0*ones(n);
    xmax[30]   = 90*ones(n);
    xmin[31]   = -0*ones(n);
    xmax[31]   = 90*ones(n);
    xmin[32]   = -0*ones(n);
    xmax[32]   = 90*ones(n);
    xmin[33]   = -0*ones(n);
    xmax[33]   = 90*ones(n);
    xmin[34]   = -0*ones(n);
    xmax[34]   = 90*ones(n);
    xmin[35]   = -0*ones(n);
    xmax[35]   = 90*ones(n);

    if 39 >= prob_k >= 36
        xmin[36]   = -1*ones(n);xmin[36][25:28] .= 0;
        xmax[36]   = +1*ones(n);
        xmin[37]   = -1*ones(n);xmin[37][25:28] .= 0;
        xmax[37]   = +1*ones(n);
        xmin[38]   = -1*ones(n);xmin[38][25:28] .= 0;
        xmax[38]   = +1*ones(n);
        xmin[39]   = -1*ones(n);xmin[39][25:28] .= 0;
        xmax[39]   = +1*ones(n);
    end
    if 46 >= prob_k >= 40
        xmin[40]   = -1*ones(n);xmin[40][27:34] .= 0;
        xmax[40]   = +1*ones(n);
        xmin[41]   = -1*ones(n);xmin[41][27:34] .= 0;
        xmax[41]   = +1*ones(n);
        xmin[42]   = -1*ones(n);xmin[42][27:34] .= 0;
        xmax[42]   = +1*ones(n);
        xmin[43]   = -1*ones(n);xmin[43][27:34] .= 0;
        xmax[43]   = +1*ones(n);
        xmin[44]   = -1*ones(n);xmin[44][27:34] .= 0;
        xmax[44]   = +1*ones(n);
        xmin[45]   = -1*ones(n);xmin[45][27:34] .= 0;
        xmax[45]   = +1*ones(n);
        xmin[46]   = -1*ones(n);xmin[46][27:34] .= 0;
        xmax[46]   = +1*ones(n);
    end
    if 49 >= prob_k >= 47
        xmin[47]   = -1*ones(n);xmin[47][11:12] .= 0;xmin[47][13:18] .= 0;
        xmax[47]   = +1*ones(n);xmax[47][11:12] .= 2;xmax[47][13:18] .= 500;
        xmin[48]   = -1*ones(n);xmin[48][11:12] .= 0;xmin[48][13:18] .= 0;
        xmax[48]   = +1*ones(n);xmax[48][11:12] .= 2;xmax[48][13:18] .= 500;
        xmin[49]   = -1*ones(n);xmin[49][11:12] .= 0;xmin[49][13:18] .= 0;
        xmax[49]   = +1*ones(n);xmax[49][11:12] .= 2;xmax[49][13:18] .= 500;
    end

    xmin[50]    = [10, 10, 35, 35, 125, 130];
    xmax[50]  = [125, 150, 210, 225, 315, 325];


    xmin_ = xmin[prob_k]
    xmax_ = xmax[prob_k]
    return n, fn, g, h, xmin_, xmax_
end

#=
function test()
    for i in 1:50
        @show i Cal_par(i)
    end
end

test()
=#

