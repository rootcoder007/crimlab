# k-core of each node

The coreness of a node is the largest k such that it belongs to a
subgraph in which every node has at least k neighbours (Seidman 1983),
computed in O(nodes + edges) by the algorithm of Batagelj and Zaversnik.
The cores are nested: the (k + 1)-core lies inside the k-core
(`proof("ntcore")`). Edges are taken as undirected, without weights.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## Usage

``` r
ntcore(from, to, nodes = NULL)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

## Value

A data frame: `node`, `core`.

## References

Seidman, S. B. (1983). Network structure and minimum degree. *Social
Networks*, 5(3), 269-287.
[doi:10.1016/0378-8733(83)90028-X](https://doi.org/10.1016/0378-8733%2883%2990028-X)

Batagelj, V. and Zaversnik, M. (2003). An O(m) algorithm for cores
decomposition of networks.
[doi:10.48550/arXiv.cs/0310049](https://doi.org/10.48550/arXiv.cs/0310049)

## Examples

``` r
ntcore(c(1, 1, 1, 2, 2, 3, 4), c(2, 3, 4, 3, 4, 4, 5))
#>   node core
#> 1    1    3
#> 2    2    3
#> 3    3    3
#> 4    4    3
#> 5    5    1
```
