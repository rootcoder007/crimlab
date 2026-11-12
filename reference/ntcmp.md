# Connected components

Weak components join nodes linked by edges in either direction; strong
components (Tarjan's algorithm, without recursion) join nodes that reach
each other along the edges' directions. Undirected, the two coincide.
Components are numbered 1, 2, ... in the order of their first node.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## Usage

``` r
ntcmp(from, to, nodes = NULL, dir = FALSE)
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

A data frame: `node`, `weak`, `strong`.

## References

Tarjan, R. (1972). Depth-first search and linear graph algorithms. *SIAM
Journal on Computing*, 1(2), 146-160.
[doi:10.1137/0201010](https://doi.org/10.1137/0201010)

## Examples

``` r
ntcmp(c("a", "b", "c", "x"), c("b", "c", "a", "y"), dir = TRUE)
#>   node weak strong
#> 1    a    1      1
#> 2    b    1      1
#> 3    c    1      1
#> 4    x    2      2
#> 5    y    2      3
```
