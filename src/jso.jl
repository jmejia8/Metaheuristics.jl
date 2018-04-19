using Distributions

function getFw(F::Float64, nfes::Int, max_nfes::Int)
    if nfes < 0.2 * max_nfes
        return 0.7F
    elseif nfes < 0.4 * max_nfes
        return 0.8F
    end
    
    return 1.2F
end

function getPbests(P::Matrix{Float64}, f::Vector{Float64}, p::Float64)
    i = sortperm(f)
    N = length(f)

    n = p*N
    if n <= 1
        n = 2
    end

    ids = i[1:round(Int, n) ]

    return ids
end

function getXr(P::Matrix{Float64}, A, i::Int)
    N  = size(P, 1)
    NA = length(A)

    r1 = rand(1:N,1)[1]
    r2 = rand(1:N+NA,1)[1]

    while r1 == i
        r1 = rand(1:N,1)[1]
    end

    while r2 == i || r1 == r2
        r2 = rand(1:N+NA,1)[1]
    end

    x1 = P[r1,:]
    x2 = 1

    if r2 > N
        r2 = 1 + ((r2-1)% N )
        x2 = A[r2]
    else
        x2 = P[r2,:]
    end

    if length(x2) == 1
        x2 = x2[1]
    end

    return x1, x2
end

function meanWL(S::Vector{Float64}, w::Vector{Float64})
    return dot(w, S.^2) / dot(w, S)
end

function jso(fobj::Function, D::Int; limits = [-100.0, 100.0])
    # conf
    p_max = 0.25
    p_min = 0.25/2
    p = p_min
    N = round(Int, 25log(D)*sqrt(D))
    N_init = N
    N_min = 4
    H = 5


    max_nfes = 10000D

    a, b = limits

    P = a + (b - a) .* rand(N, D)
    f = zeros(N)

    for i = 1:N
        f[i] = fobj(P[i,:])
    end
    nfes = N

    MF  = 0.5ones(H)
    MCR = 0.8ones(H)
    intex_M = 1

    stop = false

    g = 0

    A = []
    while !stop
        SCR   = Array{Float64}([])
        SF    = Array{Float64}([])
        Sdiff = Array{Float64}([])

        Pnext = copy(P)
        fnext = copy(f)

        N_old = N

        # println(N, " ", a, " ", b)

        pbests = getPbests(P, f, p)

        for i = 1:N
            r = rand(1:H, 1)[1]

            if r == H
                MF[r] = 0.9
                MCR[r] = 0.9
            end

            if MCR[r] < 0
                CR = 0
            else
                CR = rand(Normal(MCR[r], 0.1))
            end

            if g < 0.25max_nfes
                CR = max(CR, 0.7)
            elseif g < 0.5max_nfes
                CR = max(CR, 0.6)
            end

            F = rand(Cauchy(MF[r], 0.1))

            if g < 0.6max_nfes && F > 0.7
                F = 0.7
            end

            x = P[i,:]
            fx = f[i]


            Fw = getFw(F, nfes, max_nfes)

            pbst = rand(pbests,1)[1]
            while nfes < 0.50*max_nfes && pbst == i
                pbst = rand(pbests,1)[1]
            end

            x_pbest  = P[ pbst ,:]
            xr1, xr2 = getXr(P, A, i)
            
            v = x + Fw *(x_pbest - x) + F * (xr1 - xr2)
            u = x

            j_r = rand(1:D,1)[1]
            for j = 1:D
                if rand() <= CR || j == j_r
                    u[j] = v[j]
                end
            end


            fu = fobj(u)
            nfes += 1

            if fu <= fx
                Pnext[i,:] = u
                fnext[i] = fu
            end

            if fu < fx
                push!(A, x)
                push!(SCR, CR)
                push!(SF, F)
                push!(Sdiff, abs(fx - fu))
            end

            stop = nfes >= max_nfes
            if stop
                break
            end

        end
        

        if length(SF) > 0       
            w = Sdiff / sum(Sdiff)

            MCR[intex_M] = meanWL(SCR, w)
            MF[intex_M]  = meanWL(SF, w)

            
            intex_M += 1
            if intex_M > H
                intex_M = 1
            end
        end

        p = p_max * (1.0 - 0.5*nfes /  max_nfes)

        # update population size
        N = round(Int, N_init + ( N_min - N_init ) * nfes / max_nfes )
        if N_old != N
            Ids = sortperm(fnext)[1:N]
            P = Pnext[Ids,:]
            f = fnext[Ids]
    
            if length(A) > N
                K = randperm(length(A))[1:N]
                A = A[K]
            end

        end

        g += 1
    end

    best = indmin(f)
    return P[best, :], f[indmin(f)]


end