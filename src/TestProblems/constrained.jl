"""
    constrained1(D = 5)

Constrained test problem 1. This is a simple constrained test problem with two constraints defined as:
```math
g_1(x) = (x_1 - x_2)^2 - 6
g_2(x) = x_1 + x_2 - 2
```
and the objective function is:
```math
f(x) = (x_1 - 1)^2 + (x_2 - 1)^2
```
"""
function constrained1(D = 5)
    # Objective function
    f(x) = (sum((x .- 1).^2), [sum((x .- 1).^2) - 4, sum(sin.(x .- 1)) - 1], [(x[1] - x[2])^2])


    bounds = Array([-10.0ones(D) 10.0ones(D)]')

    x = ones(D)
    return f, bounds, [generateChild(x, f(x)) ]
end


"""
    constrained2(D = 5)

Constrained test problem 2. This is a simple constrained test problem with three constraints defined as:
```math
g_1(x) = (x_1 - x_2)^2 - 6
g_2(x) = x_1 + x_2 - 2
g_3(x) = x_1 - x_2 - 2

```
and the objective function is:
```math
f(x) = (x_1 - 1)^2 + (x_2 - 1)^2
```
"""
function constrained2(D = 5)

    f(x) = (sum((x .- 1).^2), [sum((x .- 1).^2) - 4, sum(sin.(x .- 1)) - 1, (x[1] - x[2])^2 - 6], [0.0])


    bounds = Array([-10.0ones(D) 10.0ones(D)]')

    x = ones(D)
    return f, bounds, [generateChild(x, f(x)) ]
    
end

"""
    constrained3()

Constrained test problem 3. This is a simple constrained test problem only with equality  constraints `h_i(x) = 0` defined as:
```math
h_1(x) = - y_1 + y_2 + y_3 + y_4 - 1
h_2(x) = 2x_1 - y_1 + 2y_2 - 0.5y_3 + y_5 - 1
h_3(x) = 2x_2 + 2y_1 - y_2 - 0.5y_3 + y_6 - 1
```
and the objective function is:
```math
f(x) = x_1 + 2x_2 + y_1 + y_2 + 2y_3
```
"""
function constrained3()
    f(y) = begin
        x = zeros(2)
        return x[1] + 2x[2] + y[1] + y[2] + 2y[3],
                    [0.0], # g
                    [ 
                     - y[1] + y[2] + y[3] + y[4] - 1,
                     2x[1] - y[1] + 2y[2] - 0.5y[3] + y[5] - 1,
                     2x[2] + 2y[1] - y[2] - 0.5y[3] + y[6] - 1
                    ]
    end


    bounds = Array([zeros(6) 100.0ones(6)]')
    y = zeros(6)
    y[4:6] .= 1.0
    return f, bounds, [generateChild(y, f(y)) ]
end
