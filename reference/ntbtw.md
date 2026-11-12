# Betweenness of each node

The sum, over pairs of other nodes, of the share of their shortest paths
that pass through the node (Freeman 1977), by Brandes' algorithm in
O(nodes x edges) time. Paths are counted in steps: weights are tie
strengths, not lengths, and are not taken. An undirected pair is counted
once. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## Usage

``` r
ntbtw(from, to, nodes = NULL, dir = FALSE, norm = FALSE)
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

- norm:

  `TRUE` to divide by the number of pairs of other nodes,
  `(n - 1)(n - 2)` directed or `(n - 1)(n - 2) / 2` undirected, as
  `igraph::betweenness(normalized = TRUE)` does; `NaN` for fewer than 3
  nodes.

## Value

A data frame: `node`, `btw`.

## References

Brandes, U. (2001). A faster algorithm for betweenness centrality.
*Journal of Mathematical Sociology*, 25(2), 163-177.
[doi:10.1080/0022250X.2001.9990249](https://doi.org/10.1080/0022250X.2001.9990249)

Freeman, L. C. (1977). A set of measures of centrality based on
betweenness. *Sociometry*, 40(1), 35-41.
[doi:10.2307/3033543](https://doi.org/10.2307/3033543)

## Examples

``` r
ntbtw(c(1, 2, 3), c(2, 3, 4))
#>   node btw
#> 1    1   0
#> 2    2   2
#> 3    3   2
#> 4    4   0
```
