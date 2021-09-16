function calculate_hv(points,bounds,k,nSample)

    # This function is modified from the code in
    # http://www.tik.ee.ethz.ch/sop/download/supplementary/hype/

    N,M = size(points);
    if M > 2
        # Use the estimated method for three or more objectives
        alpha = zeros(1,N); 
        for i = 1 : k 
            alpha[i] = prod((k-[1:i-1])./(N-[1:i-1]))./i; 
        end
        Fmin = ideal(points);
        # S    = unifrnd(repmat(Fmin,nSample,1),repmat(bounds,nSample,1));
        S = Fmin + (bounds - Fmin) .* rand(nSample, M)
        PdS  = zeros(Bool, N,nSample);
        dS   = zeros(1,nSample);
        for i = 1 : N
            x        = sum(points[i,:]' .- S <= 0, dims = 2) == M;
            PdS[i,x] = true;
            dS[x]    = dS[x] + 1;
        end
        F = zeros(1,N);
        for i = 1 : N
            F[i] = sum(alpha[dS[PdS[i,:]]]);
        end

        return F .* prod(bounds-Fmin) / nSample;
    end

    # Use the accurate method for two objectives
    pvec  = 1:size(points,1);
    alpha = zeros(1,k);
    for i = 1 : k 
        j = 1 : i-1; 
        alpha[i] = prod((k-j)./(N-j))./i;
    end

    return hypesub(N,points,M,bounds,pvec,alpha,k);
end

function hypesub(l,A,M,bounds,pvec,alpha,k)
    # The recursive function for the accurate method

    h     = zeros(1,l); 
    #[S,mask] = sortrows(A,M); 
    S = sortslices( view(A, :, M), dims=1 )

    # FIXME: this works but improvement is needed for performance purpurses
    mask = [ findfirs(a -> a = S[i,:]) for i in 1:size(S,1) ]

    pvec  = pvec[mask]; 
    for i = 1 : size(S,1) 
        if i < size(S,1) 
            extrusion = S[i+1,M] - S[i,M]; 
        else
            extrusion = bounds[M] - S[i,M];
        end
        if M == 1
            if i > k
                break; 
            end
            if alpha >= 0
                h[pvec[1:i]] = h[pvec[1:i]] + extrusion*alpha[i]; 
            end
        elseif extrusion > 0
            h += extrusion*hypesub(l,S[1:i,:],M-1,bounds,pvec[1:i],alpha,k); 
        end
    end

    return h
end
