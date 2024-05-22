using Documenter, Metaheuristics
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "references.bib"), style=:authoryear)

makedocs(
         bib,
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true",
                                  assets = ["assets/favicon.ico", "assets/extra_styles.css"],
                                  analytics = "UA-184071594-1",
                                  collapselevel = 1,
                                  ansicolor=true,
                                 ),
         sitename="Metaheuristics.jl",
         authors = "Jesús Mejía",
         pages = [
                  "Index" => "index.md",
                  "Tutorials" => [
                                  "tutorials/simple-tutorial.md",
                                  "tutorials/create-metaheuristic.md",
                                  "tutorials/parallelization.md",
                                  "tutorials/n-queens.md",
                                 ],
                  "Examples" =>  "examples.md",
                  "Algorithms" => [
                                  "algorithms/index.md",
                                  "algorithms/singleobjective.md",
                                  "algorithms/multiobjective.md",
                                  "algorithms/combinatorial.md",
                                 ],
                  #"Algorithms" => "algorithms.md",
                  "Problems" => "problems.md",
                  "Performance Indicators" => "indicators.md",
                  "Multi-Criteria Decision Making" => "mcdm.md",
                  "Visualization" => "visualization.md",
                  "API References" => "api.md",
                  "FAQ" => "faq.md",
                  "Contributing" => "contributing.md",
                  "References" => "references.md",
                 ]
        )



deploydocs(
           repo = "github.com/jmejia8/Metaheuristics.jl.git",
          )
