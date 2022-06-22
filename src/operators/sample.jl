function _lhs(nsamples, dim)
    # initial sample
    X = reshape([v for i in 1:dim for v in shuffle(1.0:nsamples)], nsamples, dim)
    # smooth and normalize sample
    (X - rand(nsamples, dim)) / nsamples
end

_score_lhs(M) = minimum(pairwise_distances(M))
function _scale_sample(X, bounds)
    a = view(bounds, 1,:)'
    b = view(bounds, 2,:)'
    # scale sample
    a .+ (b - a) .* X
end

function sample(method::LatinHypercubeSampling, bounds = zeros(0,0))
    X = _lhs(method.N, method.dim)
    score = _score_lhs(X)
    # sample improving
    for i in 1:method.iterations
        XX = _lhs(method.N, method.dim)
        sc = _score_lhs(XX)
        sc < score && continue
        X = XX
        score = sc
    end
    isempty(bounds) && (return X)
    _scale_sample(X, bounds)
end

#=
function sample(method::Grid, bounds = zeros(0,0))
    
    vals = Iterators.product((1:N for _ in 1:dim)...)
    X = zeros(length(vals), method.dim)
    k = 1
    for (i,j) in vals
        X[i,j] = 1
    end


    if isempty(bounds)
        a = b = range(0, 1, length=method.N)
    end

end
=#
