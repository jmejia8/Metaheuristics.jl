struct WeightedSum end
struct Tchebysheff end
struct AchievementScalarization end

eval_scalatization(F, w, params::WeightedSum) = sum(F .* w', dims=2) |> vec

function eval_scalatization(F, w, params::Tchebysheff)
    maximum(abs.(F .- minimum(F, dims=1)) .* w', dims=2) |> vec
end

function eval_scalatization(F, w, params::AchievementScalarization)
    w = copy(w)
    w[w .== zero(eltype(w))] .= eps()
    maximum((F .- minimum(F, dims=1)) ./ w', dims=2) |> vec
end


"""
    CompromiseProgramming(scalarizing)

Perform compromise programming by using the `scalarizing` function provided.
Current implemented scalarizing function are

* `WeightedSum`
* `Tchebysheff`
* `AchievementScalarization`.

### Example

```julia-repl
julia> f, bounds, pf = Metaheuristics.TestProblems.ZDT1();

julia> res = optimize(f, bounds, NSGA2());

julia> w = [0.5, 0.5];

julia> sol = best_alternative(res, w, CompromiseProgramming(Tchebysheff()))
(f = [0.38493217206706115, 0.38037042164979956], g = [0.0], h = [0.0], x = [3.849e-01, 7.731e-06, …, 2.362e-07])

julia> sol = best_alternative(res, w, CompromiseProgramming(WeightedSum()))
(f = [0.2546059308425166, 0.4958366970021401], g = [0.0], h = [0.0], x = [2.546e-01, 2.929e-06, …, 2.224e-07])

julia> sol = best_alternative(res, w, CompromiseProgramming(AchievementScalarization()))
(f = [0.38493217206706115, 0.38037042164979956], g = [0.0], h = [0.0], x = [3.849e-01, 7.731e-06, …, 2.362e-07])

julia> idx = decisionmaking(res, w, CompromiseProgramming(Tchebysheff()))
3
```
"""
mutable struct CompromiseProgramming{T, V} <: AbstractDecisionMakingMethod
    scalarizing::T
    values::Vector{V}
end

function CompromiseProgramming(scalarizing; values=zeros(0))
    CompromiseProgramming(scalarizing, values)
end


function decisionmaking(population::AbstractMatrix, w, method::CompromiseProgramming)
    values = eval_scalatization(population, w, method.scalarizing)
    method.values = values
    argmin(values)
end

function decisionmaking(population::AbstractArray{<: AbstractMultiObjectiveSolution}, args...)
    decisionmaking(fvals(population), args...)
end
