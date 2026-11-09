function initialize_closest_weight_vectors!(parameters)
    # neighborhoods given by the user
    isempty(parameters.B) || return parameters.B

    N = parameters.N
    T = parameters.T
    λ = parameters.λ

    dist = zeros(N)
    d = similar(λ[1])
    ix = Vector{Int}(undef, N)
    for i in 1:N
        λi = λ[i]
        for j in 1:N
            if j == i
                # i is never its own neighbor (even with duplicated weights)
                dist[j] = -1.0
                continue
            end
            d .= λi .- λ[j]
            dist[j] = norm(d)
        end
        # ties are broken by index, as in the stable sortperm; I[1] == i
        I = partialsortperm!(ix, dist, 1:(T + 1))
        push!(parameters.B, I[2:end])
    end
    return parameters.B
end

function update_reference_point!(z::Vector{Float64}, F::AbstractVector{<:Real})
    @inbounds for i in eachindex(z)
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
    @inbounds for i in eachindex(z)
        v = λ[i] * abs(fx[i] - z[i])
        if m < v
            m = v
        end
    end

    return m
end

# g for subproblem j: weights in column j of Λ, objectives in fx (vector) or column j of F
@inline function _g_column(fx::Vector{Float64}, Λ::Matrix{Float64}, j, z)
    m = -Inf
    @inbounds for i in eachindex(z)
        v = Λ[i, j] * abs(fx[i] - z[i])
        if m < v
            m = v
        end
    end
    return m
end

@inline function _g_column(F::Matrix{Float64}, Λ::Matrix{Float64}, j, z)
    m = -Inf
    @inbounds for i in eachindex(z)
        v = Λ[i, j] * abs(F[i, j] - z[i])
        if m < v
            m = v
        end
    end
    return m
end
