# Eigenvector centrality of each node

The leading eigenvector of the (weighted) adjacency matrix of the
undirected graph, scaled to a maximum of 1 (Bonacich 1972). It is found
by power iteration on `A + c I`, `c` the largest weight, which has the
same eigenvectors and whose leading eigenvalue is strictly dominant even
on a bipartite graph (a co-offending network between two roles, a star),
where plain power iteration oscillates. The iteration starts from 1 for
every node and stops when no value changes by `tol` or more between two
steps, or after `maxit` steps (then `conv` is `FALSE`).

## Usage

``` r
ntevc(from, to, w = NULL, nodes = NULL, tol = 1e-13, maxit = 1e+05)
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

- tol:

  The stopping tolerance: one number between 0 and 1.

- maxit:

  The most steps: one whole number, at least 1.

## Value

A list: `cent`, a data frame with `node` and `evc`; `value`, the leading
eigenvalue; `iter`, the steps taken; `conv`, whether the tolerance was
met.

## Details

On a disconnected graph the leading eigenvector lives on the components
with the largest leading eigenvalue; every other component gets 0 (each
component is iterated and scaled on its own, and a component whose
eigenvalue falls short of the largest by more than 1e-10 of it is set to
0). When several components share the largest eigenvalue each is scaled
to a maximum of 1: any mix of them is an eigenvector too, so compare
nodes within a component. A graph without edges gives 1 for every node.
Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## References

Bonacich, P. (1972). Factoring and weighting approaches to status scores
and clique identification. *Journal of Mathematical Sociology*, 2(1),
113-120.
[doi:10.1080/0022250X.1972.9989806](https://doi.org/10.1080/0022250X.1972.9989806)

## Examples

``` r
ntevc(c("hub", "hub", "hub", "a"), c("a", "b", "c", "b"))
#> $cent
#>   node       evc
#> 1  hub 1.0000000
#> 2    a 0.8546377
#> 3    b 0.8546377
#> 4    c 0.4608111
#> 
#> $value
#> [1] 2.170086
#> 
#> $iter
#> [1] 33
#> 
#> $conv
#> [1] TRUE
#> 
```
