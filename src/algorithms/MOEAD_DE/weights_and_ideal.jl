function initialize_weight_vectors!(parameters, problem)
    values = gen_weights(parameters.nobjectives, parameters.H)

    parameters.N = length(values)

    parameters.λ = values 
end

function initialize_closest_weight_vectors!(parameters, problem)
    distances = zeros(parameters.N, parameters.N)
    λ = parameters.λ
    for i in 1:parameters.N
        for j in (i+1):parameters.N
            distances[i, j] = norm(λ[i] - λ[j])
            distances[j, i] = distances[i, j]
        end
        I = sortperm(distances[i, :])
        push!(parameters.B, I[2:parameters.T+1])
    end
end

function update_reference_point!(z::Vector{Float64}, F::Vector{Float64})
    for i in 1:length(z)
        if z[i] > F[i]
            z[i] = F[i]
        end
    end
end

function update_reference_point!(z::Vector{Float64}, sol::xFgh_indiv)
    update_reference_point!(z, sol.f)
end

function update_reference_point!(z::Vector{Float64}, population)
    for sol in population
        update_reference_point!(z, sol)
    end
end

@inline g(fx, λ, z) = maximum(λ .* abs.(fx - z))



"""
    gen_weights(n_objectives, H)

"""
function gen_weights(a, b)
    nobj = a;
    H    = b;
    a    = zeros(nobj);
    d    = H;
    w    = [];
    produce_weight!(a, 1, d, H, nobj, w)
    return w
end

function  produce_weight!(a, i, d, H, nobj, w)
    for k=0:d
        if i<nobj
            a[i] = k;
            d2   = d - k;
            produce_weight!(a, i+1, d2, H, nobj, w);
        else
            a[i] = d;
            push!(w, a/H)
            break;
        end
    end
end

