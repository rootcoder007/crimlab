# Degree and strength of each node

A graph is an edge list: `from[k]` is tied to `to[k]`. The rules below
hold for every network function in crimlab.

## Usage

``` r
ntdeg(from, to, w = NULL, nodes = NULL, dir = FALSE)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- w:

  `NULL`, or one positive weight per edge (a tie's strength, such as the
  number of co-offences).

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

- dir:

  `TRUE` for a directed graph (`from` sends a tie to `to`).

## Value

A data frame with one row per node: `node`, `deg` (distinct neighbours)
and `wdeg` (summed weight); directed, `node`, `in`, `out`, `deg`
(`in + out`), `win`, `wout` and `wdeg`. The degrees sum to twice the
number of distinct edges: see `proof("ntdeg")`.

## Details

- Nodes are identified by their labels, compared as text
  ([`as.character()`](https://rdrr.io/r/base/character.html); a factor
  by its labels, unused levels ignored). Results have one row per node,
  in the order nodes first appear: each edge's `from`, then its `to`,
  edge by edge. With `nodes`, the order is that of `nodes`, which must
  hold every label of the edges once; labels in `nodes` that no edge
  names are isolated nodes.

- A self-loop (`from[k] == to[k]`) is dropped; its node stays.

- Repeated edges (the same pair, in either direction when undirected)
  are merged into one edge whose weight is the sum of their weights (1
  each when `w` is `NULL`). So `deg` counts distinct neighbours and
  `wdeg` counts edges or their weight.

- A missing value in any argument is an error.

## References

Wasserman, S. and Faust, K. (1994). *Social Network Analysis: Methods
and Applications*. Cambridge University Press.
[doi:10.1017/CBO9780511815478](https://doi.org/10.1017/CBO9780511815478)

## Examples

``` r
ntdeg(c("a", "a", "b", "a"), c("b", "c", "c", "b"))
#>   node deg wdeg
#> 1    a   2    3
#> 2    b   2    3
#> 3    c   2    2
ntdeg(c("a", "a"), c("b", "c"), nodes = c("c", "b", "a", "z"), dir = TRUE)
#>   node in out deg win wout wdeg
#> 1    c  1   0   1   1    0    1
#> 2    b  1   0   1   1    0    1
#> 3    a  0   2   2   0    2    2
#> 4    z  0   0   0   0    0    0
```
