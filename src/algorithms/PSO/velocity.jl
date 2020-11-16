@inline function velocity(x, v, pbest, gbest, parameters)
    r1 = parameters.C1 * rand()
    r2 = parameters.C2 * rand()
    parameters.ω * v + r1 * (pbest - x) + r2 * (gbest - x)
end
