# N-Queens

In the N-queen problem, $N$ queens must be placed on an $N\times N$ chessboard without interfering with each other.
There can never be more than one queen on a chess row, column or diagonal because of the queen's ability to move vertically, horizontally and diagonally. 

## Solution Representation

A naive representation is to use a binary matrix to represent a chessboard $c_{ij}$
where $c_{ij}=1$ would indicate a queen in row $i$ and column $j$; however the search
space cardinality (number of solutions) is huge for this representation.
A permutation-based representation can be used instead.
A permutation will contain in position $j$ the row $i$ where a queen is, removing a lot of
infeasible solutions (queens attacking horizontally and vertically).

Let's compute the number of solutions (a.k.a. cardinality):

```@repl
using Metaheuristics
N = 8 # for 8-queens
cardinality(BitArrays(N*N))
cardinality(Permutations(N))
```

## Objective Function

It is necessary to minimize the number of attacks, the following objective function is
proposed to this end:

```@example queens
function attacks(chessboard)
    N = length(chessboard)
    n_attacks = 0
    # check attack in both diagonas for each queen
    for i = 1:N
        Δrows = i + chessboard[i]
        Δcols = i - chessboard[i]
        for j = (i+1):N
            # check diagonal [\]
            n_attacks +=  j + chessboard[j] == Δrows ? 1 : 0
            # check inverse diagonal [/]
            n_attacks +=  j - chessboard[j] == Δcols ? 1 : 0
        end
    end
    2n_attacks
end
```

It can observed that horizontal and vertical are not considered due to the
adopted solution representation.


## Let's find a solution

To minimize the number of attacks, let's use a Genetic Algorithm:

```@example queens
using Metaheuristics # hide
N = 8
optimize(attacks, Permutations(N), GA); # hide
result = optimize(attacks, Permutations(N), GA)
```

It can be observed that the number of attacks is:

```@example queens
n_attacks = minimum(result)
```

The optimal permutation is:

```@example queens
perm = minimizer(result)
```

Corresponding chessboard configuration:

```@example queens
chessboard = zeros(Bool, N, N)
for (i,j) = enumerate(perm); chessboard[i,j] = true;end
chessboard
```

## The 20-queens case

```@example queens
N = 20
perm = optimize(attacks, Permutations(N), GA) |> minimizer
```

```@example queens
attacks(perm)
```


```@example queens
chessboard = zeros(Bool, N, N)
for (i,j) = enumerate(perm); chessboard[i,j] = true;end
chessboard
```
