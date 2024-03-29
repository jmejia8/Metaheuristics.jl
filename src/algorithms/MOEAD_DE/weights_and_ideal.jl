function initialize_closest_weight_vectors!(parameters)
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

function g(fx, λ, z) 
    m = -Inf
    for i in eachindex(z)
        v = λ[i] * abs(fx[i] - z[i])
        if m < v
            m = v
        end
    end
    
    return m
end

