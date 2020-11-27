using Documenter, Metaheuristics

makedocs(
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets = ["assets/favicon.ico"],
                                  analytics = "G-MCJ59XDJNN",
                                 ),
         sitename="Metaheuristics.jl",
         authors = "Jesús Mejía",
         pages = [
                  "Introduction" => "index.md",
                  "Tutorial" => "tutorial.md",
                  "Algorithms" => "algorithms.md",
                  "Problems" => "problems.md",
                  "Performance Idicators" => "indicators.md",
                  "Visualization" => "visualization.md",
                  "API References" => "api.md",
                  "Contributing" => "contributing.md",
                 ]
        )



deploydocs(
           repo = "github.com/jmejia8/Metaheuristics.jl.git",
          )
