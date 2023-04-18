function sample(
        sampler::Sampler{M, S},
        n = 100
    ) where {M <: RandomSampler, S <: AtomicSearchSpace}

    X = rand(sampler.method.rng, sampler.searchspace, n)
    # convert to matrix
    [X[i][j] for i in eachindex(X), j in eachindex(first(X))]
end


function sample(
        sampler::Sampler{M, S}, args...
    ) where {M <: GridSampler, S <: AtomicSearchSpace}

    X = collect(sampler)
    # convert to matrix
    [X[i][j] for i in eachindex(X), j in eachindex(first(X))]
end

function sample(
        ::Type{T},
        search_space::AtomicSearchSpace,
        n = 100,
    ) where T <: AbstractSampler

    sample(T(search_space), n)
end

function sample(search_space::AtomicSearchSpace, n = 100)

    if cardinality(search_space) <= n
        sampler = GridSampler(search_space)
    else
        sampler = RandomSampler(search_space)
    end

    sample(sampler, n)
end


