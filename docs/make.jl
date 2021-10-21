using Documenter, Metaheuristics
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"))

makedocs(
         bib,
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets = ["assets/favicon.ico"],
                                  analytics = "UA-184071594-1",
                                 ),
         sitename="Metaheuristics.jl",
         authors = "Jesús Mejía",
         pages = [
                  "Index" => "index.md",
                  "Tutorials" => ["tutorials/simple-tutorial.md", "tutorials/create-metaheuristic.md"],
                  "Examples" =>  "examples.md",
                  "Algorithms" => "algorithms.md",
                  "Problems" => "problems.md",
                  "Performance Indicators" => "indicators.md",
                  "Visualization" => "visualization.md",
                  "API References" => "api.md",
                  "Contributing" => "contributing.md",
                  "References" => "references.md",
                 ]
        )



deploydocs(
           repo = "github.com/jmejia8/Metaheuristics.jl.git",
          )
