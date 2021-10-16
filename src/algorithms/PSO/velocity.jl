function PSO_velocity(x, v, pbest, gbest, C1, C2, ω)
    r1 = C1 * rand()
    r2 = C2 * rand()
    ω * v + r1 * (pbest - x) + r2 * (gbest - x)
end

function velocity(x, v, pbest, gbest, parameters)
    PSO_velocity(x, v, pbest, gbest, parameters.C1, parameters.C2, parameters.ω)
end

