function get_last_front(population, N = length(population))
    let f::Int = 1
        ind = 0
        _indnext = findlast(x -> x.rank == f, population)
        indnext = isnothing(_indnext) ? 0 : _indnext

        while 0 < indnext <= N
            ind = indnext
            f += 1
            indnext = findlast(x -> x.rank == f, population)
        end

        indnext == 0 && (indnext = length(population))

        return ind+1:indnext
    end
    
end

function update_front_contribution!(population, N = length(population))
    # FIXME performance improvement required (update front)
    fast_non_dominated_sort!(population)

    last_front = get_last_front(population)

    compute_contribution!(view(population, last_front))
end


function compute_contribution!(population)
    if isempty(population)
        # nothing to do
        return 
    end

    M = length(fval(population[1]))
    N = length(population)

    ΔS = fill(Inf, N) 

    # objetive function values
    Fs = fvals(population)

    if M == 2
        rank = sortslices(Fs, dims=1)

        for i = 2 : N-1
            ΔS[rank[i]] = (Fs[rank[i+1],1] - Fs[rank[i],1]) .* (Fs[rank[i-1],2] - Fs[rank[i],2]);
        end

    elseif N > 1
        ΔS = calculate_hv(Fs, nadir(population)*1.1, 1, 10000)
    end

    worst = argmin(ΔS)
    deleteat!(population, worst)
    
end

