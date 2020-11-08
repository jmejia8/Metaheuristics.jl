using Documenter, Metaheuristics

makedocs(
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets = ["assets/favicon.ico"]
                                 ),
         sitename="Metaheuristics.jl",
         authors = "Jesús Mejía",
         pages = [
                  "Introduction" => "index.md",
                  "Tutorial" => "tutorial.md",
                  "Algorithms" => "algorithms.md",
                  "API References" => "api.md",
                 ]
        )



deploydocs(
           repo = "github.com/jmejia8/Metaheuristics.jl.git",
          )
