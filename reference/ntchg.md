# Cheeger bounds on a network's bottleneck

The conductance of a set of nodes S is the weight of the edges leaving
it over the smaller of its volume and that of the rest (volume: summed
strength); the graph's conductance h is the least over all sets, a
measure of how cleanly the network splits in two. Finding h is hard, but
the second-smallest eigenvalue `lambda2` of the normalised Laplacian
`I - D^-1/2 A D^-1/2` brackets it: `lambda2 / 2 <= h <= sqrt(2 lambda2)`
(Cheeger's inequality; the lower half is `proof("ntchg")`). The sweep
over the eigenvector finds a set whose conductance (`sweep`) is at most
`sqrt(2 lambda2)`; h is at most `sweep`.

## Usage

``` r
ntchg(from, to, w = NULL, nodes = NULL, tol = 1e-10, maxit = 1e+05)
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

A list: `bound`, a named vector `lambda2`, `lower` (`lambda2 / 2`),
`upper` (`sqrt(2 lambda2)`) and `sweep`; `set`, the labels of the sweep
set; `iter`; `conv`.

## Details

`lambda2` is found by power iteration on `(I + D^-1/2 A D^-1/2) / 2`
kept orthogonal to `D^1/2 1`, from a fixed pseudo-random start, until
the residual is below `tol` or after `maxit` steps (`conv` `FALSE`).
Convergence is slow when `lambda2` and the next eigenvalue are close
(long paths, large sparse graphs). The graph is undirected and must be
connected; weights are tie strengths. Graph rules as in
[`ntdeg()`](https://rootcoder007.github.io/crimlab/reference/ntdeg.md).

## References

Chung, F. R. K. (1997). *Spectral Graph Theory*. CBMS Regional
Conference Series in Mathematics 92. American Mathematical Society.
[doi:10.1090/cbms/092](https://doi.org/10.1090/cbms/092)

## Examples

``` r
# two triangles joined by one edge
ntchg(c(1, 1, 2, 4, 4, 5, 3), c(2, 3, 3, 5, 6, 6, 4))
#> $bound
#>   lambda2     lower     upper     sweep 
#> 0.2046664 0.1023332 0.6397911 0.1428571 
#> 
#> $set
#> [1] "4" "5" "6"
#> 
#> $iter
#> [1] 31
#> 
#> $conv
#> [1] TRUE
#> 
```
