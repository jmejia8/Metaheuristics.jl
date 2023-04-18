function PSO_velocity(x, v, pbest, gbest, C1, C2, ω, rng)
    r1 = C1 * rand(rng)
    r2 = C2 * rand(rng)
    ω * v + r1 * (pbest - x) + r2 * (gbest - x)
end

function velocity(x, v, pbest, gbest, parameters, rng)
    PSO_velocity(x, v, pbest, gbest, parameters.C1, parameters.C2, parameters.ω, rng)
end

