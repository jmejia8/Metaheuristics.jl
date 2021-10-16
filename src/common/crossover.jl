function gen_β(β, η, D, R)
    α = 2.0 .- β .^ (-  η - 1.0 )
    mask = R .<= 1.0 ./ α
    s = 1.0 / (η + 1.0)
    βq = [ mask[i] ?  (R[i] * α[i])^s : (1.0 / (2.0 - R[i]*α[i]))^s for i in 1:D]
    βq
end

function SBX_crossover(vector1, vector2, bounds, η=15, p_variable = 0.9)
    xu = view(bounds, 2,:)
    xl = view(bounds, 1,:)
    D = length(vector1)

    do_crossover = ones(Bool, D)
    do_crossover[rand(D) .> p_variable] .= false
    do_crossover[ abs.( vector2 - vector1 ) .<= eps() ] .= false

    y1 = min.( vector1, vector2 )
    y2 = max.( vector1, vector2 )
    Δ = max.(eps(), y2 - y1)


    R = rand(D)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c1 = @. 0.5*(y1 + y2 -  βq*Δ)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D, R) 
    c2 = @. 0.5*(y1 + y2 +  βq*Δ)

    # swap
    mask = rand(Bool, D)
    cc = copy(c1)
    c1[mask] = c2[mask]
    c2[mask] = cc[mask]

    cc1 = copy(vector1)
    cc1[do_crossover] = c1[do_crossover]
    cc2 = copy(vector2)
    cc2[do_crossover] = c2[do_crossover]


    reset_to_violated_bounds!(cc1, bounds)
    reset_to_violated_bounds!(cc2, bounds)

    return cc1, cc2
end

function DE_crossover(x, u, CR)
    D = length(x)
    # binomial crossover
    v = zeros(D)
    j_rand = rand(1:D)

    # binomial crossover
    for j = 1:D
        if rand() < CR || j == j_rand
            v[j] = u[j]
        else
            v[j] = x[j]
        end
    end

    return v
end

