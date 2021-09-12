function  CEC2021_func_1_25(x,func)
    # CEC2021_func Real world Multi-objective Constrained Optimization Test Suite 
    # Abhishek Kumar (email: abhishek.kumar.eee13@iitbhu.ac.in, Indian Institute of Technology (BHU), Varanasi) 

    # x -----> ps X D where 'ps': number of population and 'D': Dimension of
    # the problem.
    # f -----> Objective Function Value.
    # g -----> Inequality Consstraints Value; ps X ng where 'ng': number of
    # inequality constraints.
    # h -----> Equality Constraints Value; ps X nh where 'nh': number of
    # equality constraints.
    # prob_k -> Index of problem.

    D = size(x);

    ## Meachnical problems

    if func == 1
        ## Pressure Vessal Problems
        x1 = round(x[1]);
        x2 = round(x[2]);
        x3 = x[3];
        x4 = x[4];
        z1 = 0.0625*x1;
        z2 = 0.0625*x2;
        ## objective function
        f = zeros(2);
        f[1] = 1.7781 .* z1 .* x3 .^ 2+0.6224 .* z1 .* x2 .* x4+3.1661 .* z1 .^ 2 .* x4+19.84 .* z1 .^ 2 .* x3;
        f[2] = -pi .* x3 .^ 2 .* x4-(4/3) .* pi .* x3 .^ 3;
        ## constraints
        g = zeros(2);
        h = zeros(1);
        g[1] = 0.00954 .* x3-z2;
        g[2] = 0.0193 .* x3-z1;
    elseif func == 2
        ## Vibrating Platform
        d1 = x[1];
        d2 = x[2];
        d3 = x[3];
        b = x[4];
        L = x[5];
        rho1 = 100; rho2 = 2770; rho3 = 7780;
        E1 = 1.6; E2 = 70; E3 = 200;
        c1 = 500; c2 = 1500; c3 = 800;
        mu = 2*b .* (rho1 .* d1+rho2 .* (d2-d1)+rho3 .* (d3-d2));
        EI = (2*b ./ 3) .* (E1 .* d1 .^ 3+E2 .* (d2 .^ 3-d1 .^ 3)+rho3 .* (d3-d2));
        ## objective function
        f = zeros(2);
        f[1] = (-pi) ./ (2*L) .^ 2 .* (abs(EI ./ mu)) .^ 0.5;
        f[2] = 2*b .* L .* (c1 .* d1+c2 .* (d2-d1)+c3 .* (d3-d2));
        ## constraints
        g = zeros(5);
        h = zeros(1);
        g[1] = mu .* L -2800;
        g[2] = d1-d2;
        g[3] = d2-d1-0.15;
        g[4] = d2-d3;
        g[5] = d3-d2-0.01;
    elseif func == 3
        ## Two bar Truss design Problems
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        ## objective function
        f = zeros(2);
        f[1] = x1 .* (16+x3 .^ 2) .^ (0.5)+x2 .* (1+x3 .^ 2) .^ (0.5);
        f[2] = (20 .* (16+x3 .^ 2) .^ (0.5)) ./ (x3 .* x1);
        ## constraints
        g = zeros(3);
        h = zeros(1);
        g[1] = f[1]-0.1;
        g[2] = f[2]-1e5;
        g[3] = (80 .* (1+x3 .^ 2) .^ (0.5)) ./ (x3 .* x2)-1e5;
    elseif func == 4
        ## Weldan Beam Design Problem
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        P = 6000;
        L = 14;
        E = 30e6;
        tmax = 13600;
        sigmax = 30000;
        G = 12e6;
        Pc = (4.013 .* E .* ((x3 .^ 2+x4 .^ 6) ./ 36) .^ 0.5) ./ (L .^ 2) .* (1-x3 ./ (2 .* L) .* (E ./ (4*G))^(0.5));
        sigma = (6 .* P .* L) ./ (x4 .* x3 .^ 2);
        J = 2 .* (sqrt(2) .* x1 .* x2 .* (x2 .^ 2/12+((x1+x3) ./ 2) .^ 2));
        R = sqrt(x2 .^ 2/4+((x1+x3) ./ 2) .^ 2);
        M = P .* (L+x2 ./ 2);
        tho1 = P ./ (sqrt(2) .* x1 .* x2);
        tho2 = M .* R ./ J;
        tho = sqrt(tho1 .^ 2+2 .* tho1 .* tho2 .* x2 ./ (2*R)+tho2 .^ 2);
        ## objective function
        f = zeros(2);
        f[1] = 1.10471 .* x1 .^ 2 .* x2+0.04811 .* x3 .* x4 .* (14+x2);
        f[2] = (4 .* P .* L .^ 3) ./ (E .* x4 .* x3 .^ 3);
        ## constraints
        g = zeros(4);
        h = zeros(1);
        g[1] = tho -tmax;
        g[2] = sigma - sigmax;
        g[3] = x1-x4;
        g[4] = P - Pc;
    elseif func == 5
        ## Disc Brake Design Problem
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        ## objective function 
        f = zeros(2);
        f[1] = 4.9e-5 .* (x2 .^ 2-x1 .^ 2) .* (x4-1);
        f[2] = 9.82e6 .* (x2 .^ 2-x1 .^ 2) ./ (x3 .* x4 .* (x2 .^ 3-x1 .^ 3));
        ## constraints
        g = zeros(4);
        h = zeros(1);
        g[1] = 20 - (x2 - x1);
        g[2] = x3 ./ (3.14 .* (x2 .^ 2-x1 .^ 2)) - 0.4;
        g[3] = 2.22e-3 .* x3 .* (x2 .^ 3-x1 .^ 3) ./ (x2 .^ 2-x1 .^ 2) .^ 2-1;
        g[4] = 900 - 2.66e-2 .* x3 .* x4 .* (x2 .^ 3-x1 .^ 3) ./ (x2 .^ 2-x1 .^ 2);
    elseif func == 6
        ## Speed Reducer Design Problem
        x1 = x[1];
        x2 = x[2];
        x3 = round(x[3]);
        x4 = x[4];
        x5 = x[5];
        x6 = x[6];
        x7 = x[7];
        ## objective function
        f = zeros(2);
        f[1] = 0.7854 .* x1 .* x2 .^ 2 .* (10 .* x3 .^ 2/3+14.933 .* x3-43.0934) - 1.508 .* x1 .* (x6 .^ 2+x7 .^ 2)+7.477 .* (x6 .^ 3+x7 .^ 3)+0.7854 .* (x4 .* x6 .^ 2+x5 .* x7 .^ 2);
        f[2] = sqrt((745 .* x4 ./ (x2 .* x3)) .^ 2+1.69e7) ./ (0.1 .* x6 .^ 3);
        ## constraints
        g = zeros(11);
        h = zeros(1);
        g[1] = 1 ./ (x1 .* x2 .^ 2 .* x3)-1/27;
        g[2] = 1 ./ (x1 .* x2 .^ 2 .* x3 .^ 2)-1/397.5;
        g[3] = x4 .^ 3 ./ (x2 .* x3 .* x6 .^ 4)-1/1.93;
        g[4] = x5 .^ 3 ./ (x2 .* x3 .* x7 .^ 4)-1/1.93;
        g[5] = x2 .* x3-40;
        g[6] = x1 ./ x2-12;
        g[7] = -x1 ./ x2+5;
        g[8] = 1.9-x4+1.5 .* x6;
        g[9] = 1.9-x5+1.1 .* x7;
        g[10] = f[2]-1300;
        g[11] = sqrt((745 .* x5 ./ (x2 .* x3)) .^ 2+1.575e8) ./ (0.1*x7 .^ 3)-850;
    elseif func == 7
        ## Gear Train Design Problems
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        ## objective function
        f = zeros(2);
        f[1] = abs(6.931-x3 .* x4 ./ (x1 .* x2));
        f[2] = maximum(x);
        ## constraints
        g = zeros(1);
        h = zeros(1);
        g[1] = f[1] ./ 6.931-0.5;
    elseif func == 8
        ## Car Side Impact Desifn Problem
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        x5 = x[5];
        x6 = x[6];
        x7 = x[7];
        VMBP = 10.58-0.674 .* x1 .* x2-0.67275 .* x2;
        VFD  = 16.45-0.489 .* x3 .* x7-0.843 .* x5 .* x6;
        ## objective function
        f = zeros(3);
        f[1] = 1.98+4.9 .* x1 .* 6.67 .* x2+6.98 .* x3+4.01 .* x4+1.78 .* x5+1e-5 .* x6+2.73 .* x7;
        f[2] = 4.72-0.5 .* x4-0.19 .* x2 .* x3;
        f[3] = 0.5 .* (VMBP+VFD);
        ## constraints
        g = zeros(9);
        h = zeros(1);
        g[1] = -1+1.16-0.3717 .* x2 .* x4-0.0092928 .* x3;
        g[2] = -0.32+0.261-0.0159 .* x1 .* x2-0.06486 .* x1-0.019 .* x2 .* x7+0.0144 .* x2 .* x5+0.0154464 .* x6;
        g[3] = -0.32+0.74-0.61 .* x2-0.031296 .* x3-0.031872 .* x7+0.227 .* x2 .^ 2;
        g[4] = -0.32+0.214+0.00817 .* x5-0.045195 .* x1-0.0135168 .* x1+0.03099 .* x2 .* x6-0.018 .* x2 .* x7+0.007176 .* x3+0.023232 .* x3-0.00364 .* x5 .* x6-0.018 .* x2 .^ 2;
        g[5] = -32+33.86+2.95 .* x3-5.057 .* x1 .* x2-3.795 .* x2-3.4431 .* x7+1.45728;
        g[6] = -32+28.98+3.818 .* x3-4.2 .* x1 .* x2+1.27296 .* x6-2.68065 .* x7;
        g[7] = -32+46.36-9.9 .* x2-4.4505 .* x1;
        g[8] = f[2]-4;
        g[9] = VMBP - 9.9;
    elseif func == 9
        ## Four Bar Plane Truss
        F = 10; E = 2e5; L = 200; sig = 10;
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        ## objective function
        f = zeros(2);
        f[1] = L .* (2*x1+sqrt(2) .* x2+sqrt(2) .* x3+x4);
        f[2] = F .* L ./ E .* (2 ./ x1+2 .* sqrt(2) ./ x2-2 .* sqrt(2) ./ x3+2 ./ x4);
        ## constraints
        g = zeros(1);
        h = zeros(1);
    elseif func == 10
        ## Two bar plane Truss
        x1 = x[1];
        x2 = x[2];
        rho = 0.283;
        h = 100;
        P = 104;
        E = 3e7;
        rho0 = 2e4;
        ## objective function
        f = zeros(2);
        f[1] = 2 .* rho .* h .* x2 .* sqrt(1+x1 .^ 2);
        f[2] = rho .* h .* (1+x1 .^ 2) .^ 1.5 .* (1+x1 .^ 4) .^ 0.5 ./ (2 .* sqrt(2) .* E .* x1 .^ 2 .* x2);
        ## constraints
        g = zeros(2);
        h = zeros(2);
        g[1] = P .* (1+x1) .* (1+x1 .^ 2) .^ 0.5 ./ (2 .* sqrt(2) .* x1 .* x2)-rho0;
        g[2] = P .* (-x1+1) .* (1+x1 .^ 2) .^ 0.5 ./ (2 .* sqrt(2) .* x1 .* x2)-rho0;
    elseif func == 11
        ## Water Resource Management Problem
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        ## objectives
        f = zeros(5);
        f[1] = 106780.37  .*  (x2 + x3) + 61704.67 ;
        f[2] = 3000  .*  x1 ;
        f[3] = 305700  .*  2289  .*  x2  ./  ((0.06 .* 2289) ^ 0.65) ;
        f[4] = 250  .*  2289  .*  exp(-39.75 .* x2+9.9 .* x3+2.74) ;
        f[5] = 25  .*  (1.39  ./ (x1 .* x2) + 4940 .* x3 -80) ;
        ## Constraints   
        g = zeros(7);
        h = zeros(1);
        g[1] = 1 - (0.00139 ./ (x1 .* x2)+4.94 .* x3-0.08);
        g[2] = 1 - (0.000306 ./ (x1 .* x2)+1.082 .* x3-0.0986);
        g[3] = 50000 - (12.307 ./ (x1 .* x2) + 49408.24 .* x3+4051.02);
        g[4] = 16000 - (2.098 ./ (x1 .* x2)+8046.33 .* x3-696.71);
        g[5] = 10000 - (2.138 ./ (x1 .* x2)+7883.39 .* x3-705.04);
        g[6] = 2000 - (0.417 .* x1 .* x2 + 1721.26 .* x3-136.54);
        g[7] = 550 - (0.164 ./ (x1 .* x2)+631.13 .* x3-54.48);
        g = -g;
    elseif func == 12
        ## Simply supported I-beam design  
        x1 = x[1];
        x2 = x[2];
        x3 = x[3];
        x4 = x[4];
        P = 600;
        L = 200;
        E = 2e4;
        ## objectives
        f = zeros(2);
        f[1] = 2  .*  x2  .*  x4 + x3  .*  (x1- 2 .* x4);
        f[2] = P .* L .^ 3 ./ (48 .* E .* (x3  .* ((x1 - 2 .* x4) .^ 3)+2 .* x2 .* x4 .* (4 .* x4 .* x4+3 .* x1 .* (x1-2 .* x4))) ./ 12);
        ## constraints
        g = zeros(1);
        h = zeros(1);
        g[1] = -16 + 180000*x1 ./ (x3 .* ((x1 - 2*x4) .^ 3)+2*x2 .* x4 .* (4 .* x4 .* x4+3*x1 .* (x1-2*x4)))
        + 15000 .* x2 ./ ((x1-2 .* x4) .* x3 .^ 3+2 .* x4 .* x2 .^ 3);
    elseif func == 13
        ## Gear box design 
        x1 = x[1];
        x2 = x[2];
        x3 = round(x[3]);
        x4 = x[4];
        x5 = x[5];
        x6 = x[6];
        x7 = x[7];
        ## objective function
        f = zeros(3);
        f[1] = 0.7854 .* x1 .* x2 .^ 2 .* (10 .* x3 .^ 2/3+14.933 .* x3-43.0934) - 1.508 .* x1 .* (x6 .^ 2+x7 .^ 2)+7.477 .* (x6 .^ 3+x7 .^ 3)+0.7854 .* (x4 .* x6 .^ 2+x5 .* x7 .^ 2);
        f[2] = sqrt((745 .* x4 ./ (x2 .* x3)) .^ 2+1.69e7) ./ (0.1 .* x6 .^ 3);
        f[3] = sqrt((745 .* x5 ./ (x2 .* x3)) .^ 2+1.575e8) ./ (0.1*x7 .^ 3);
        ## constraints
        g = zeros(11);
        h = zeros(1);
        g[1] = 1 ./ (x1 .* x2 .^ 2 .* x3)-1/27;
        g[2] = 1 ./ (x1 .* x2 .^ 2 .* x3 .^ 2)-1/397.5;
        g[3] = x4 .^ 3 ./ (x2 .* x3 .* x6 .^ 4)-1/1.93;
        g[4] = x5 .^ 3 ./ (x2 .* x3 .* x7 .^ 4)-1/1.93;
        g[5] = x2 .* x3-40;
        g[6] = x1 ./ x2-12;
        g[7] = -x1 ./ x2+5;
        g[8] = 1.9-x4+1.5 .* x6;
        g[9] = 1.9-x5+1.1 .* x7;
        g[10] = f[2]-1300;
        g[11] = f[3]-1100;
    elseif func == 14
        ## Multiple-disk clutch brake design problem 
        Mf = 3; Ms = 40; Iz = 55; n = 250; Tmax = 15; s = 1.5; delta = 0.5; 
        Vsrmax = 10; rho = 0.0000078; pmax = 1; mu = 0.6; Lmax = 30; delR = 20;
        Rsr = 2 ./ 3 .* (x[2] .^ 3-x[1] .^ 3) ./ (x[2] .^ 2 .* x[1] .^ 2);
        Vsr = pi .* Rsr .* n ./ 30;
        A   = pi .* (x[2] .^ 2-x[1] .^ 2);
        Prz = x[4] ./ A;
        w   = pi .* n ./ 30;
        Mh  = 2/3 .* mu .* x[4] .* x[5] .* (x[2] .^ 3-x[1] .^ 3) ./ (x[2] .^ 2-x[1] .^ 2);
        T   = Iz .* w ./ (Mh+Mf);
        ## objective function
        f = zeros(2);
        f[1] = pi .* (x[2] .^ 2-x[1] .^ 2) .* x[3] .* (x[5]+1) .* rho;
        f[2] = T;
        ## constraints
        g = zeros(8)
        g[1] = -x[2]+x[1]+delR;
        g[2] = (x[5]+1) .* (x[3]+delta)-Lmax;
        g[3] = Prz-pmax;
        g[4] = Prz .* Vsr-pmax .* Vsrmax;
        g[5] = Vsr-Vsrmax;
        g[6] = T-Tmax;
        g[7] = s .* Ms-Mh;
        g[8] = -T;
        h      = zeros(1);
    elseif func == 15
        ## Spring design Problem  
        x1 = round(x[1]);
        x2 = x[2];
        d  = [0.009,0.0095,0.0104,0.0118,0.0128,0.0132,0.014,
              0.015, 0.0162, 0.0173, 0.018, 0.020, 0.023, 0.025,
              0.028, 0.032, 0.035, 0.041, 0.047, 0.054, 0.063,
              0.072, 0.080, 0.092, 0.0105, 0.120, 0.135, 0.148,
              0.162, 0.177, 0.192, 0.207, 0.225, 0.244, 0.263,
              0.283, 0.307, 0.331, 0.362,0.394,0.4375,0.500];
        x3 = d[max(1,min(42,round(Int, x[3])))];# x3 = x3(:);
        ## constants
        cf = (4 .* x2 ./ x3-1) ./ (4 .* x2 ./ x3-4)+0.615 .* x3 ./ x2;
        K  = (11.5 .* 10 .^ 6 .* x3 .^ 4) ./ (8 .* x1 .* x2 .^ 3);
        lf = 1000 ./ K + 1.05 .* (x1+2) .* x3;
        sigp = 300 ./ K;
        ## objective function
        f = zeros(2);
        f[1] = (pi .^ 2 .* x2 .* x3 .^ 2 .* (x1+2)) ./ 4;
        f[2] = (8000 .* cf .* x2) ./ (pi .* x3 .^ 3);
        ## constraints
        g = zeros(8)
        g[1] = (8000 .* cf .* x2) ./ (pi .* x3 .^ 3)-189000;
        g[2] = lf-14;
        g[3] = 0.2-x3;
        g[4] = x2-3;
        g[5] = 3-x2 ./ x3;
        g[6] = sigp - 6;
        g[7] = sigp+700 ./ K+1.05 .* (x1+2) .* x3-lf;
        g[8] = 1.25-700 ./ K;
        h = zeros(1);

    elseif func == 16
        ## Cantilever beam design problem 
        x1 = x[1];
        x2 = x[2];

        P = 1;
        E = 207000000;
        Sy = 300000;
        delta_max = 0.005;
        rho = 7800;

        ## objectives
        f = zeros(2);
        f[1] = 0.25  .*  rho  .*  pi  .*  x2  .*  x1 .^ 2;
        f[2] = (64  .*  P  .*  x2 .^ 3) ./ (3  .*  E  .*  pi  .*  x1 .^ 4);

        ## constraints
        g = zeros(2);
        h = zeros(1);
        g[1] = -Sy + (32  .*  P  .*  x2) ./ (pi  .*  x1 .^ 3);
        g[2] = -delta_max + (64  .*  P  .*  x2 .^ 3) ./ (3  .*  E  .*  pi  .*  x1 .^ 4);

    elseif func == 17
        ## Bulk carriers design problem 
        L = x[1];
        B = x[2];
        D = x[3];
        T = x[4];
        V_k = x[5];
        C_B = x[6];

        a = 4977.06 .* C_B .^ 2 - 8105.61 .* C_B + 4456.51;
        b = -10847.2 .* C_B .^ 2 + 12817 .* C_B - 6960.32;
        F_n = 0.5144 ./ (9.8065  .*  L) .^ 0.5;
        P = ((1.025 .* L .* B .* T .* C_B) .^ (2/3) .* V_k .^ 3) ./ (a + b .* F_n);

        W_s = 0.034 .* L .^ 1.7 .* B .^ 0.6 .* D .^ 0.4 .* C_B .^ 0.5;
        W_o = L .^ 0.8 .* B .^ 0.6 .* D .^ 0.3 .* C_B .^ 0.1;
        W_m = 0.17 .* P .^ 0.9;
        ls = W_s+W_o+W_m;

        D_wt = 1.025 .* L .* B .* T .* C_B-ls;
        F_c = 0.19 .* 24 .* P ./ 1000 + 0.2;
        D_cwt = D_wt - F_c .* ((5000 .* V_k) ./ 24 + 5)-2 .* D_wt .^ 0.5;
        R_trp = 350 ./ ((5000 .* V_k) ./ 24 + 2 .* (D_cwt ./ 8000 + 0.5));
        ac = D_cwt .* R_trp;
        S_d = 5000 .* V_k ./ 24;

        C_c = 0.2 .* 1.3  .*  (2000 .* W_s .^ 0.85 + 3500 .* W_o + 2400 .* P .^ 0.8);
        C_r = 40000 .* D_wt .^ 0.3;
        C_v = (1.05 .* 100 .* F_c .* S_d + 6.3 .* D_wt .^ 0.8) .* R_trp;

        ## objectives
        f = zeros(3);
        f[1]= (C_c + C_r + C_v) ./ ac;
        f[2] = ls;
        f[3] = -ac;
        ## constraints
        g = zeros(9);
        h = zeros(1);
        g[1] = L ./ B - 6;
        g[2]= 15 - L ./ D;
        g[3] = 19 - L ./ T;
        g[4] = 0.45 .* D_wt .^ 0.31 - T;
        g[5] = 0.7 .* D + 0.7 - T;
        g[6] = 0.32 - F_n;
        g[7] = 0.53 .* T + ((0.085 .* C_B - 0.002) .* B .^ 2) ./ (T .* C_B)-(1 + 0.52 .* D) - 0.07 .* B;
        g[8] = D_wt - 3000;
        g[9] = 500000 - D_wt;
        g = -g;
    elseif func == 18
        ## Front rail design problem 
        hh = x[1];
        w = x[2];
        t = x[3];

        Ea = 14496.5;
        Fa = 234.9;
        E = -70973.4 + 958.656  .*  w+ 614.173  .*  hh - 3.827  .*  w  .*  hh + 57.023  .*  w  .*  t + 63.274  .*  hh  .*  t-3.582  .*  w .^ 2 - 1.4842  .*  hh .^ 2 - 1890.174  .*  t .^ 2;
        F = 111.854 - 20.210  .*  w + 7.560  .*  hh - 0.025  .*  w  .*  hh + 2.731  .*  w  .*  t - 1.479  .*  hh  .*  t + 0.165  .*  w .^ 2;

        ## objectives
        f = zeros(2);
        f[1] = Ea ./ E;
        f[2] = F ./ Fa;

        ## constraints
        g = zeros(3);
        h = zeros(1);
        g[1] = (hh - 136) .*  (146 - hh);
        g[2] = (w - 58)  .*  (66 - w);
        g[3] = (t - 1.4)  .*  (2.2 - t);
        g = -g;
    elseif func == 19
        ## Multi-product batch plant
        ## constant
        S = [2 3 4;
             4 6 3];
        t = [ 8 20 8;
             16 4 4];
        H = 6000; alp = 250; beta = 0.6;
        Q1 = 40000; Q2 = 20000;
        ## decision Variable
        N1 = round(x[1]); N2 = round(x[2]); N3 = round(x[3]);
        V1 = x[4]; V2 = x[5]; V3 = x[6];
        TL1 = x[7]; TL2 = x[8];
        B1 = x[9]; B2 = x[10];
        ## objective function
        f = zeros(3);
        f[1] = alp .* (N1 .* V1 .^ beta+N2 .* V2 .^ beta+N3 .* V3 .^ beta);
        f[2] = 65 .* (Q1 ./ B1+Q2 ./ B2)+0.08 .* Q1+0.1 .* Q2;
        f[3] = Q1 .* TL1 ./ B1+Q2 .* TL2 ./ B2;
        ## constraints
        g = zeros(10)
        g[1] = Q1 .* TL1 ./ B1+Q2 .* TL2 ./ B2-H;
        g[2] = S[1,1] .* B1+S[2,1] .* B2-V1;
        g[3] = S[1,2] .* B1+S[2,2] .* B2-V2;
        g[4] = S[1,3] .* B1+S[2,3] .* B2-V3;
        g[5] = t[1,1]-N1 .* TL1;
        g[6] = t[1,2]-N2 .* TL1;
        g[7] = t[1,3]-N3 .* TL1;
        g[8] = t[2,1]-N1 .* TL2;
        g[9] = t[2,2]-N2 .* TL2;
        g[10] = t[2,3]-N3 .* TL2;
        h = zeros(1);
    elseif func == 20
        ## Hydro-static thrust bearing design problem
        R = x[1]; Ro = x[2];  mu = x[3]; Q = x[4];
        gamma = 0.0307; C = 0.5; n = -3.55; C1 = 10.04;
        Ws = 101000; Pmax = 1000; delTmax = 50; hmin = 0.001;
        gg = 386.4; N = 750;
        P    = (log10(log10(8.122*1e6 .* mu+0.8))-C1) ./ n;
        delT = 2 .* (10 .^ P-560);
        Ef   = 9336 .* Q .* gamma .* C .* delT;
        h    = (2 .* pi .* N ./ 60) .^ 2 .* 2 .* pi .* mu ./ Ef .* (R .^ 4 ./ 4-Ro .^ 4 ./ 4)-1e-5;
        Po   = (6 .* mu .* Q ./ (pi .* h .^ 3)) .* log(R ./ Ro);
        W    = pi .* Po ./ 2 .* (R .^ 2-Ro .^ 2) ./ (log(R ./ Ro)-1e-5);
        ##  objective function
        f = zeros(2);
        f[1] = (Q .* Po ./ 0.7+Ef) ./ 12;
        f[2] = gamma ./ (gg .* Po) .* (Q ./ (2 .* pi .* R .* h));
        ##  constraints
        g = zeros(7)
        g[1] = Ws-W;
        g[2] = Po-Pmax;
        g[3] = delT-delTmax;
        g[4] = hmin-h;
        g[5] = Ro-R;
        g[6] = f[2]-0.001;
        g[7] = W ./ (pi .* (R .^ 2-Ro .^ 2)+1e-5)-5000;
        h      = zeros(1);
    elseif func == 21
        ## Crash energy management for high-speed train
        x1 = x[1]; x2 = x[2]; x3 = x[3];
        x4 = x[4]; x5 = x[5]; x6 = x[6];
        ## objective function
        f = zeros(2);
        f[1] = 1.3667145844797-0.00904459793976106*x1-0.0016193573938033*x2-0.00758531275221425*x3-0.00440727360327102*x4-0.00572216860791644*x5-0.00936039926190721*x6+2.62510221107328*10^(-6)*(x1 .^ 2)+4.92982681358861*10^(-7)*(x2 .^ 2)+2.25524989067108*10^(-6)*(x3 .^ 2)+
        1.84605439400301*10^(-6)*(x4 .^ 2)+2.17175358243416*10^(-6)*(x5 .^ 2)+3.90158043948054*10^(-6)*(x6 .^ 2)+4.55276994245781*10^(-7)*x1 .* x2-6.37013576290982*10^(-7)*x1 .* x3+8.26736480446359*10^(-7)*x1 .* x4+5.66352809442276*10^(-8)*x1 .* x5-3.20213897443278*10^(-7)*x1 .* x6+
        1.18015467772812*10^(-8)*x2 .* x3+9.25820391546515*10^(-8)*x2 .* x4-1.05705364119837*10^(-7)*x2 .* x5-4.74797783014687*10^(-7)*x2 .* x6-5.02319867013788*10^(-7)*x3 .* x4+9.54284258085225*10^(-7)*x3 .* x5+1.80533309229454*10^(-7)*x3 .* x6-1.07938022118477*10^(-6)*x4 .* x5-
        1.81370642220182*10^(-7)*x4 .* x6-2.24238851688047*10^(-7)*x5 .* x6; 
        f[2] = -1.19896668942683+3.04107017009774*x1+1.23535701600191*x2+2.13882039381528*x3+2.33495178382303*x4+2.68632494801975*x5+3.43918953617606*x6-7.89144544980703*10^(-4)*(x1 .^ 2)-2.06085185698215*10^(-4)*(x2 .^ 2)-7.15269900037858*10^(-4)*(x3 .^ 2)-7.8449237573837*10^(-4)*(x4 .^ 2)-
        9.31396896237177*10^(-4)*(x5 .^ 2)-1.40826531972195*10^(-3)*(x6 .^ 2)-1.60434988248392*10^(-4)*x1 .* x2+2.0824655419411*10^(-4)*x1 .* x3-3.0530659653553*10^(-4)*x1 .* x4-8.10145973591615*10^(-5)*x1 .* x5+6.94728759651311*10^(-5)*x1 .* x6+1.18015467772812*10^(-8)*x2 .* x3+
        9.25820391546515*10^(-8)*x2 .* x4-1.05705364119837*10^(-7)*x2 .* x5+1.69935290196781*10^(-4)*x2 .* x6+2.32421829190088*10^(-5)*x3 .* x4-2.0808624041163476*10^(-4)*x3 .* x5+1.75576341867273*10^(-5)*x3 .* x6+2.68422081654044*10^(-4)*x4 .* x5+4.39852066801981*10^(-5)*x4 .* x6+
        2.96785446021357*10^(-5)*x5 .* x6;
        ## constraints
        g = zeros(4);
        h = zeros(1);
        g[1] = f[1]-5;
        g[2] = -f[1];
        g[3] = f[2] - 28;
        g[4] = -f[2];
        ## chemical problems
    elseif func == 22
        ## Haverly's Pooling Problem
        x1 = x[1]; x2 = x[2]; x3 = x[3];
        x4 = x[4]; x5 = x[5]; x6 = x[6];
        x7 = x[7]; x8 = x[8]; x9 = x[9];
        ## objective function
        f = zeros(2);
        f[1] = -9*x1-15*x2+6*x3+16*x4;
        f[2] = 10 .* (x5+x6);
        ## constraints
        g = zeros(2);
        h = zeros(4);
        g[1] = x9 .* x7+2*x5-2.5*x1;
        g[2] = x9 .* x8+2*x6-1.5*x2;
        h[1] = x7+x8-x4-x3;
        h[2] = x1-x5-x7;
        h[3] = x2-x6-x8;
        h[4] = x9 .* x7+x9 .* x8-3 .* x3-x4;
    elseif func == 23
        ## Reactor Network Design
        k1 = 0.09755988; k2 = 0.99*k1;
        k3 = 0.0391908; k4 = 0.9*k3;
        x1 = x[1]; x2 = x[2]; x3 = x[3];
        x4 = x[4]; x5 = x[5]; x6 = x[6];
        ## objective function
        f = zeros(2);
        f[1] = -x4;
        f[2] = x5 .^ (0.5)+x6 .^ (0.5);
        ## constraints
        g = zeros(1);
        h = zeros(4);
        g[1] = f[2]-4;
        h[1] = k1 .* x5 .* x2 + x1 -1;
        h[2] = k3 .* x5 .* x3+x3+x1-1;
        h[3] = k2 .* x6 .* x2 - x1 + x2;
        h[4] = k4 .* x6 .* x4 + x2-x1+x4-x3;
    elseif func == 24
        ## Heat Exchanger Network Design
        x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
        x6 = x[6]; x7 = x[7]; x8 = x[8]; x9 = x[9];
        ## objective function
        f = zeros(3);
        f[1] = 35 .* x1 .^ (0.6)+ 35 .* x2 .^ (0.6);
        f[2] = 200 .* x1 .* x4-x3;
        f[3] = 200 .* x2 .* x6-x5;
        ## constraints
        g = zeros(1);
        h = zeros(6);
        h[1] = x3 - 1e4 .* (x7-100);
        h[2] = x5 - 1e4 .* (300-x7);
        h[3] = x3 - 1e4 .* (600-x8);
        h[4] = x5 - 1e4 .* (900-x9);
        h[5] = x4 .* log(abs(x8-100)+1e-6)-x4 .* log(abs(600-x7)+1e-6)-x8+x7+500;
        h[6] = x6 .* log(abs(x9-x7)+1e-6)-x6 .* log(600)-x9+x7+600;
        ## process Design and Synthesis
    elseif func == 25
        ## process synthesis problem
        x1 = x[1]; x2 = round(x[2]);
        ## objective function
        f = zeros(2);
        f[1] = x2 + 2*x1;
        f[2] = -x1 .^ 2 -x2;
        ## constraints
        g = zeros(2);
        h = zeros(1);
        g[1] = f[2] +1.25;
        g[2] = x1 + x2 -1.6;

    end

    return f, g, h
end
