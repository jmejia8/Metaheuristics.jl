include("calc-hv.jl")

function update_population!(population, offspring, n_samples)
    push!(population, offspring)

    # FIXME performance improvement required (update front)
    fast_non_dominated_sort!(population) 

    last_front = get_last_front(population)

    compute_contribution!(population, last_front, n_samples)
end

function compute_contribution!(population, last_front, n_samples)
    if isempty(population)
        # nothing to do
        return 
    end

    M = length(fval(population[1]))


    # objetive function values
    Fs = fvals(population[last_front])
    N = size(Fs, 1)
    ΔS = fill(Inf, N) 

    if M == 2
        rank = sortslicesperm(Fs, dims=1)

        for i = 2 : N-1
            ΔS[rank[i]] = (Fs[rank[i+1],1] - Fs[rank[i],1]) .* (Fs[rank[i-1],2] - Fs[rank[i],2]);
        end

    elseif N > 1
        ΔS = calculate_hv(Fs, nadir(population)*1.1, 1, n_samples)
    end

    # save contribution for further analysis
    for i in eachindex(ΔS)
        population[last_front[i]].crowding = ΔS[i]
    end


    worst = argmin(ΔS)
    deleteat!(population, last_front[worst])
end

get_contribution(sol::xFgh_indiv) = sol.crowding

function get_last_front(population, N = length(population) - 1)
    let f::Int = 1
        ind = 0
        _indnext = findlast(x -> x.rank == f, population)
        indnext = isnothing(_indnext) ? 0 : _indnext

        while 0 < indnext <= N
            ind = indnext
            f += 1
            _indnext = findlast(x -> x.rank == f, population)
            indnext = isnothing(_indnext) ? 0 : _indnext
        end

        indnext == 0 && (indnext = length(population))

        return ind+1:indnext
    end
    
end

