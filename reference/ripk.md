# Ripley's K function

Ripley's K function of a point pattern in the rectangle `win` at radii
`r`, and Besag's `L = sqrt(K / pi)`. Under complete spatial randomness
`K(r) = pi r^2` and `L(r) = r`. Four estimators: `"none"` (no edge
correction), `"border"` (only points at least `r` from the edge are
centres), `"modif"` (the modified border estimator, with the area of the
eroded window; `NaN` once `r` reaches half a side) and `"iso"` (Ripley's
isotropic correction, each pair weighted by the inverse of the share of
the circle through the other point that lies in the rectangle, computed
exactly). The normalisation is `n (n - 1)`, as in spatstat.

## Usage

``` r
ripk(x, y = NULL, win, r, corr = "iso", metric = "euclid")
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y` (an unnamed two-column
  matrix is read as `x` then `y`). At least two points, all inside
  `win`; no missing values.

- win:

  The study rectangle `c(xmin, xmax, ymin, ymax)`, in the unit of the
  coordinates. The intensity is the number of points over its area.

- r:

  The radii: non-negative numbers in strictly increasing order, in the
  unit of the coordinates.

- corr:

  The edge correction: `"iso"`, `"border"`, `"modif"` or `"none"`.

- metric:

  The distance: `"euclid"` or `"manhattan"` (`|dx| + |dy|`, a street
  grid aligned with the axes), as in
  [`spwt()`](https://rootcoder007.github.io/crimlab/reference/spwt.md)
  and
  [`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md).
  Under Manhattan distance `K(r) = 2 r^2` and `L = sqrt(K / 2)` under
  randomness; the border and modified border corrections carry over (a
  Manhattan ball of radius `r` lies in the rectangle exactly when its
  centre is at least `r` from every edge), while the isotropic
  correction, the share of a circle, has no Manhattan form and is
  refused.

## Value

A data frame with one row per radius: `r`, `k` (in the square of the
coordinates' unit) and `l` (`sqrt(k / pi)`, or `sqrt(k / 2)` under
Manhattan distance).

## References

Ripley, B. D. (1976). The second-order analysis of stationary point
processes. *Journal of Applied Probability*, 13(2), 255-266.
[doi:10.2307/3212829](https://doi.org/10.2307/3212829)

Ripley, B. D. (1977). Modelling spatial patterns. *Journal of the Royal
Statistical Society, Series B*, 39(2), 172-192.
[doi:10.1111/j.2517-6161.1977.tb01615.x](https://doi.org/10.1111/j.2517-6161.1977.tb01615.x)

## Examples

``` r
set.seed(1)
ripk(runif(100), runif(100), win = c(0, 1, 0, 1), r = seq(0, 0.25, 0.05))
#>      r           k          l
#> 1 0.00 0.000000000 0.00000000
#> 2 0.05 0.007686885 0.04946525
#> 3 0.10 0.033207387 0.10281167
#> 4 0.15 0.080766297 0.16033936
#> 5 0.20 0.131964756 0.20495289
#> 6 0.25 0.201545381 0.25328618
```
