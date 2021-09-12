include("Cal_par.jl")
include("problems-1-25.jl")
include("problems-26-40.jl")
include("problems-41-50.jl")

function  CEC2021_func(x,func)
    if 1 <= func <= 25
        return CEC2021_func_1_25(x, func)
    elseif 26 <= func <= 40
        return CEC2021_func_26_40(x, func)
    elseif 41 <= func <= 50
        return CEC2021_func_41_50(x, func)
    end

end


function test()
    for i in 1:50 
        n, fn, g, h, xmin, xmax = Cal_par(i)

        x = xmin + 0.9*(xmax - xmin) 
        @show i  CEC2021_func(x, i)
    end
end

test()
