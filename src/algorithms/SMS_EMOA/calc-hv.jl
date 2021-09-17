function sortslicesperm(A::AbstractArray; dims::Union{Integer, Tuple{Vararg{Integer}}}, kws...)
    _sortslicesperm(A, Val{dims}(); kws...)
end

function _sortslicesperm(A::AbstractArray, d::Val{dims}; kws...) where dims
    itspace = Base.compute_itspace(A, d)
    vecs = map(its->view(A, its...), itspace)
    p = sortperm(vecs; kws...)
    if ndims(A) == 2 && isa(dims, Integer) && isa(A, Array)
        # At the moment, the performance of the generic version is subpar
        # (about 5x slower). Hardcode a fast-path until we're able to
        # optimize this.
        # return dims == 1 ? A[p, :] : A[:, p]
        return p
    else
        B = similar(A)
        for (x, its) in zip(p, itspace)
            B[its...] = vecs[x]
        end
        B
    end
end


function calculate_hv(points,bounds,k,nSample)

    # This function is modified from the code in
    # http://www.tik.ee.ethz.ch/sop/download/supplementary/hype/

    N,M = size(points);
    if M > 2
        # Use the estimated method for three or more objectives
        alpha = zeros(1,N); 
        for i = 1 : k 
            J = 1:i-1
            alpha[i] = prod((k .- J)./(N .- J ))./i; 
        end
        Fmin = ideal(points);
        # S    = unifrnd(repmat(Fmin,nSample,1),repmat(bounds,nSample,1));
        S = Fmin' .+ (bounds - Fmin)' .* rand(nSample, M)
        PdS  = zeros(Bool, N,nSample);
        dS   = zeros(Int, nSample);
        for i = 1 : N
            x        = sum(points[i,:]' .- S .<= 0, dims = 2) .== M;
            mask = x[:,1]
            PdS[i, mask] .= true;
            dS[mask]    = dS[mask] .+ 1;
        end
        F = zeros(N);
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
    _A = view(A, :, M)
    maks = sortslicesperm( _A, dims=1 )
    S = A[mask,:]

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
