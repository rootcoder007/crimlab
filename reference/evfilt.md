# Eigenvector spatial filtering

Adds Moran eigenvectors of the neighbour weights to the regression of
`v` on an intercept and `x` by forward selection, so that the residuals
lose their spatial autocorrelation (Griffith 2003; Tiefelsdorf and
Griffith 2007). The Moran eigenvectors are the eigenvectors of the
symmetrised weights `(W + W') / 2` orthogonal to the constant (computed
in the Helmert basis), each with Moran's I `n / S0 lambda`; the
candidates are those with Moran's I above `-1 / (n - 1)` (`"positive"`,
the usual choice), below it, or all. Selection by `method`: `"moran"`
takes the vector that most reduces the absolute residual Moran's I until
it is below `tol`; `"aic"` and `"press"` (the leave-one-out prediction
error) while they fall; `"r2"` while R-squared rises by `tol` or more.
At least two residual degrees of freedom are kept.

## Usage

``` r
evfilt(
  v,
  w,
  x = NULL,
  method = "moran",
  cand = "positive",
  tol = 0.1,
  maxv = NULL
)
```

## Arguments

- v:

  The outcome: finite numbers, one per area.

- w:

  The neighbours, as in
  [`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
  for as many areas as `v`.

- x:

  `NULL`, or covariates: a numeric matrix (or vector) with one row per
  area. An intercept is added.

- method:

  `"moran"`, `"aic"`, `"press"` or `"r2"`.

- cand:

  `"positive"`, `"negative"` or `"all"`.

- tol:

  The stopping tolerance: one positive number.

- maxv:

  `NULL`, or the largest number of vectors.

## Value

A list of class `evfilt`: `sel` (the selected eigenvectors in the order
they were added, numbered by decreasing eigenvalue), `coef` (the
intercept, `x`, then the vectors), `imoran` (the residual Moran's I),
`aic`, `r2`, `fr2` (the gain in R-squared from the filter), `resid` and
`mi` (the Moran's I of every eigenvector).

## References

Griffith, D. A. (2003). *Spatial Autocorrelation and Spatial Filtering*.
Springer, Berlin.
[doi:10.1007/978-3-540-24806-4](https://doi.org/10.1007/978-3-540-24806-4)

Tiefelsdorf, M. and Griffith, D. A. (2007). Semiparametric filtering of
spatial autocorrelation: the eigenvector approach. *Environment and
Planning A*, 39(5), 1193-1221.
[doi:10.1068/a37378](https://doi.org/10.1068/a37378)

## See also

Other spatial epidemiology:
[`bayout()`](https://rootcoder007.github.io/crimlab/reference/bayout.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
w <- 1 * (abs(outer(1:8, 1:8, "-")) == 1)
v <- c(1, 2, 2.5, 4, 3.5, 3, 1.5, 1)
evfilt(v, w, method = "aic")
#> evfilt: 1 eigenvectors selected, residual Moran's I -0.7032, R-squared 0.9199
```
