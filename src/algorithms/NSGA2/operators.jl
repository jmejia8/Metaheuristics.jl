function tournament_selection(P, i, is_better)
    a = i
    b = rand(1:length(P))
    
    is_better(P[a], P[b]) || (!is_better(P[b], P[a]) && P[a].crowding > P[b].crowding ) ? P[a] : P[b]
end


function gen_β(β, η, D)
    α = 2.0 .- β .^ (-  η - 1.0 )
    R = rand(D)
    mask = R .<= 1.0 ./ α
    s = 1.0 / (η + 1.0)
    βq = [ mask[i] ?  (R[i] * α[i])^s : (1.0 / (2 - R[i]*α[i]))^s for i in 1:D]
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


    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D) 
    c1 = @. 0.5*(y1 + y2 -  βq*Δ)

    β = @. 1.0 + (2.0 * (y1 - xl) / Δ)
    βq = gen_β(β, η, D) 
    c2 = @. 0.5*(y1 + y2 -  βq*Δ)

    # swap
    mask = rand(Bool, D)
    cc = copy(c1)
    c1[mask] = c2[mask]
    c2[mask] = cc[mask]

    cc1 = copy(vector1)
    cc1[do_crossover] = c1[do_crossover]
    cc2 = copy(vector2)
    cc2[do_crossover] = c2[do_crossover]

    return cc1, cc2
end

function polynomial_mutation!(vector, bounds, η=20, prob = 1 / length(vector))
    do_mutation = rand(length(vector)) .< prob

    xu = view(bounds, 2,do_mutation)
    xl = view(bounds, 1,do_mutation)
    x = view(vector, do_mutation)

    δ1 = (x - xl) ./ (xu - xl)
    δ2 = (xu - x) ./ (xu - xl)

    D = length(xu)
    R = rand(D)
    mask = rand(Bool, D) .< 0.5
    s = η+1.0
    mut_pow = 1.0 / (η + 1.0)
    δq = [ mask[i] ?
            ^(2R[i] + (1 - 2R[i]) * ^(1 - δ1[i], s), mut_pow) - 1.0 :
            (2.0 * (1.0 - R[i]) + 2.0 * (R[i] - 0.5) * ^(1.0 - δ2[i], s))^mut_pow
            for i in 1:D
        ]

    vector[do_mutation] = x + δq .* ( xu - xl)
    # correct using reset to bound

end
