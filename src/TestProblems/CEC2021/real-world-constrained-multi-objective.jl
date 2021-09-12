using MAT
using LinearAlgebra
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

#=
function test()
    data = matread("data0.5.mat")["st"]
    fs = data["f"]
    gs = data["g"]
    hs = data["h"]
    xs = data["x"]

    for i in 1:50 
        fval = fs[i][1,:]

        if gs[i] isa Real
            gval = [gs[i]]
        elseif isempty(gs[i])
            gval = zeros(1)
        else
            gval = gs[i][1,:]
        end


        if hs[i] isa Real
            hval = [hs[i]]
        elseif isempty(hs[i])
            hval = zeros(1)
        else
            hval = hs[i][1,:]
        end

        n, fn, g, h, xmin, xmax = Cal_par(i)

        #x = xmin + 0.5*(xmax - xmin) 
        x = xs[i][1,:]
        fx, gx, hx = CEC2021_func(x, i)


        print("fun = ", i, "  ")

        if length(fval) != length(fx)
            println(" error.")
            @show length(fval), length(fx), fn
        end

        err = norm(fval - fx)
        err1 = norm(gval - gx)
        err2 = norm(hval - hx)
        if err < 1e-10 && err1 < 1e-10 && err2 < 1e-10
            println("   OK")
        else

            @show x
            @show err
            println("")
            println("fm: ", fval)
            println("fj: ", fx)
            println("----")


            @show err1
            println("gm: ", gval)
            println("gj: ", gx)
            println("----")

            @show err2
            println("hm: ", hval)
            println("hj: ", hx)
            println("<<<<<<<<<<<<<<<<<<")
            println("<<<<<<<<<<<<<<<<<<")
        end
    end
end

test()
=#
