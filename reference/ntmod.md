# Modularity of a partition

`Q = sum over communities of in / 2m - res (tot / 2m)^2`, where `in` is
the weight of edges inside the community (counted from both ends), `tot`
the summed strength of its nodes and `2m` the total strength (Newman and
Girvan 2004), on the undirected, weighted graph. At `res = 1` it lies
between -1/2 and 1 (`proof("ntmod")`). `NaN` for a graph without edges.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## Usage

``` r
ntmod(from, to, memb, w = NULL, nodes = NULL, res = 1)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- memb:

  The community of each node: labels, either named by node (every node
  named once; other names are ignored) or unnamed, one per node in the
  order the network functions report nodes.

- w:

  `NULL`, or one positive weight per edge (a tie's strength, such as the
  number of co-offences).

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

- res:

  The resolution: one number, at least 0.

## Value

The modularity, one number.

## References

Newman, M. E. J. and Girvan, M. (2004). Finding and evaluating community
structure in networks. *Physical Review E*, 69(2), 026113.
[doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113)

## Examples

``` r
ntmod(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4),
      memb = c(1, 1, 1, 2, 2, 2))
#> [1] 0.3571429
```
