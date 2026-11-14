# Quadrat count test

Counts the points in an `nx` by `ny` grid of equal tiles of the
rectangle `win` and tests complete spatial randomness by Pearson's
chi-square, `X2 = sum (n_k - e)^2 / e` with `e` the mean count, on \`nx
ny

- 1`degrees of freedom: the upper tail against clustering and the two-sided p-value (twice the smaller tail,`quadrat.test()`'s default), which also flags regular patterns. A point on an inner tile boundary goes to the tile above or to its right, one on the outer edge to the tile inside. These are the values of `spatstat.explore::quadrat.test()\`.
  The counts feed
  [`morsta()`](https://rootcoder007.github.io/crimlab/reference/morsta.md).

## Usage

``` r
qdrtst(x, y = NULL, win, nx, ny = nx)
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y`. At least two points; no
  missing values.

- win:

  The study rectangle `c(xmin, xmax, ymin, ymax)`, holding every point.

- nx, ny:

  The tiles across and up: whole numbers from 1, at least 2 tiles.

## Value

A list of class `qdrtst`: `counts`, a matrix with one row per row of
tiles from the top and one column per column from the left; `stat`, a
named vector: `x2`, `df`, `p` (upper tail), `ptwo` and `lambda` (the
points per unit area).

## References

Greig-Smith, P. (1952). The use of random and contiguous quadrats in the
study of the structure of plant communities. *Annals of Botany*, 16(2),
293-316.
[doi:10.1093/oxfordjournals.aob.a083317](https://doi.org/10.1093/oxfordjournals.aob.a083317)

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(100, 0, 2), y = runif(100))
q <- qdrtst(p, win = c(0, 2, 0, 1), nx = 4, ny = 2)
q
#> qdrtst: X2 5.6 on 7 df, p 0.587 (two-sided 0.826)
#>      [,1] [,2] [,3] [,4]
#> [1,]   13   16   10   15
#> [2,]    7   16   11   12
morsta(q$counts)
#>       imor       mclu       muni       imst        chi          p 
#>  0.9858586  1.0910380  0.9463623 -0.1318235  5.6000000  0.5871510 
```
