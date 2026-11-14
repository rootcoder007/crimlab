# Closeness of each node

Distances are counted in steps (weights are not taken), along edges'
directions when `dir = TRUE` (out-distances). A disconnected graph is
handled in the two standard ways, both returned: `cls` and `ncls` count
only the nodes a node reaches (as `igraph::closeness()` does), so a node
close to the few it reaches scores high; `hrm`, the harmonic centrality,
sums `1 / distance` over every other node, an unreachable one adding 0,
so reaching more nodes always scores more. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## Usage

``` r
ntcls(from, to, nodes = NULL, dir = FALSE)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

- dir:

  `TRUE` for a directed graph (`from` sends a tie to `to`).

## Value

A data frame: `node`, `reach` (nodes reached), `cls` (1 / the summed
distance to them), `ncls` (`reach` / that sum, the inverse mean
distance; `cls` and `ncls` are `NaN` when nothing is reached) and `hrm`
(the harmonic centrality).

## References

Freeman, L. C. (1978). Centrality in social networks: conceptual
clarification. *Social Networks*, 1(3), 215-239.
[doi:10.1016/0378-8733(78)90021-7](https://doi.org/10.1016/0378-8733%2878%2990021-7)

Marchiori, M. and Latora, V. (2000). Harmony in the small-world.
*Physica A*, 285(3-4), 539-546.
[doi:10.1016/S0378-4371(00)00311-3](https://doi.org/10.1016/S0378-4371%2800%2900311-3)

## Examples

``` r
ntcls(c(1, 2, 3), c(2, 3, 4), nodes = 1:5)
#>   node reach       cls ncls      hrm
#> 1    1     3 0.1666667 0.50 1.833333
#> 2    2     3 0.2500000 0.75 2.500000
#> 3    3     3 0.2500000 0.75 2.500000
#> 4    4     3 0.1666667 0.50 1.833333
#> 5    5     0       NaN  NaN 0.000000
```
