include("vnd.jl")


Base.@kwdef struct VNS{I, S, L, C, N} <: AbstractVNS
    initial::I
    neighborhood_shaking::S
    neighborhood_local::L
    local_search::C
    neighborhood_change::N
end

# change neighborhood
struct SequentialChange end
struct CyclicChange end

function VNS(;initial=nothing,neighborhood_shaking=nothing, neighborhood_local=nothing,
        local_search=nothing, neighborhood_change=SequentialChange(),
        options=MH.Options(), information=MH.Information())

    parameters = VNS(initial, neighborhood_shaking, neighborhood_local,
                     local_search, neighborhood_change)

    MH.Algorithm(parameters; options, information)
end

function shake(x, neighborhood, rng)
    # select at random
    neighborhood_structure(x, neighborhood, rng)
end

function neighborhood_change(old, new, k, ::SequentialChange)
    if MH.is_better(new, old)
        return 1, new
    end
    k + 1, old
end

function neighborhood_change(old, new, k, ::CyclicChange)
    k += 1
    if MH.is_better(new, old)
        return k, new
    end
    k, old
end

function MH.update_state!(status, parameters::VNS, problem, information, options, args...; kargs...)
    # current solution
    sol = first(status.population)
    k = 1
    while k <= length(parameters.neighborhood_shaking)
        x = MH.get_position(sol)
        neighborhood = parameters.neighborhood_shaking[k]
        # select x at random from kth neighborhood
        xp = shake(x, neighborhood, options.rng)

        # perform local search around xp using VND
        # TODO: update this for considering other VNS variants (for the local search)
        vnd = VND(;initial=xp, neighborhood=parameters.neighborhood_local,
                  local_search = FirstImprovingSearch())
        _res_local = MH.optimize(problem.f, problem.search_space, vnd)
        sol_new = _res_local.best_sol

        # neighborhood change or not?
        k, sol = neighborhood_change(sol, sol_new, k, parameters.neighborhood_change)

        # save best result so far (internal use only, VNS doesn't use it)
        if MH.is_better(sol_new, status.best_sol)
            status.best_sol = sol_new
        end
        problem.f_calls += _res_local.f_calls

    end
    # save sol (x) for the next iteration of VNS
    status.population = [sol]
end

