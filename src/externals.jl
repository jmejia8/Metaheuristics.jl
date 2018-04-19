if Pkg.installed("Distributions") == nothing
    Pkg.install("Distributions")
end

using Distributions