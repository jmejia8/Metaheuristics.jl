if Pkg.installed("Distributions") == nothing
    Pkg.add("Distributions")
end

using Distributions