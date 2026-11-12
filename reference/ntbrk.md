# Brokerage roles

Counts, for each node v, the two-paths `i -> v -> j` with `i != j` and
no tie `i -> j` (v brokers between them), by the roles of Gould and
Fernandez (1989), given each node's group: `coord` (coordinator: i, v
and j in one group), `itin` (itinerant: i and j in one group, v in
another), `repr` (representative: v in i's group, j outside), `gate`
(gatekeeper: v in j's group, i outside) and `liais` (liaison: three
groups). Undirected, each tie runs both ways, so each unordered pair of
neighbours counts twice and `repr` equals `gate`. Weights are not taken.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).
Raw counts only: the expected counts and z-scores of Gould and Fernandez
are not given.

## Usage

``` r
ntbrk(from, to, grp, nodes = NULL, dir = FALSE)
```

## Arguments

- from, to:

  Labels of the two ends of each edge: vectors of the same length (zero
  when `nodes` gives every node) with no missing values.

- grp:

  The group of each node: labels, named by node or one per node in node
  order, as `memb` in
  [`ntmod()`](https://rootcoder007.github.io/crimlab/reference/ntmod.md).

- nodes:

  `NULL`, or every node's label once: the node order of the results and
  the way to include isolated nodes.

- dir:

  `TRUE` for a directed graph (`from` sends a tie to `to`).

## Value

A data frame: `node`, `grp` (as text), `coord`, `itin`, `repr`, `gate`,
`liais` and `total`.

## References

Gould, R. V. and Fernandez, R. M. (1989). Structures of mediation: a
formal approach to brokerage in transaction networks. *Sociological
Methodology*, 19, 89-126.
[doi:10.2307/270949](https://doi.org/10.2307/270949)

## Examples

``` r
ntbrk(c("a", "b", "b"), c("b", "c", "d"),
      grp = c(a = "x", b = "y", c = "y", d = "z"), dir = TRUE)
#>   node grp coord itin repr gate liais total
#> 1    a   x     0    0    0    0     0     0
#> 2    b   y     0    0    0    1     1     2
#> 3    c   y     0    0    0    0     0     0
#> 4    d   z     0    0    0    0     0     0
```
