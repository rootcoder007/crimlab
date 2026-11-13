# Knox test of near repeats

Asks whether incidents close together in space also happen close
together in time more often than chance gives: the near-repeat pattern
of burglary and shootings (Knox 1964; Townsley, Homel and Chaseling
2003). Every pair of incidents within the largest distance is put in a
band of distance, `(ds[k - 1], ds[k]]`, and a band of time,
`(dt[l - 1], dt[l]]` (the first bands start at 0 and include it), and
each cell's count is compared with its expectation when the times are
exchangeable among the locations, which is exact: `E = S_k T_l / P`,
with `S_k` the pairs in distance band `k`, `T_l` the pairs in time band
`l` and `P = n (n - 1) / 2`. Each cell has a one-sided permutation
p-value, the times shuffled among the locations `nsim` times from
`seed`, and the Poisson tail `P(X >= obs)` with the expectation as its
mean, the classical Knox approximation for a rare cell. Cells are tested
separately: a table of many cells shows some small p-values by chance
alone. A near-repeat pattern also arises when one offender commits a
series; the test speaks to timing, not to its cause.

## Usage

``` r
knxnr(x, y = NULL, t, ds, dt, nsim = 999, seed = NULL, lonlat = FALSE)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as
  [`sf::st_coordinates()`](https://r-spatial.github.io/sf/reference/st_coordinates.html)
  names them, or `lon` and `lat` with `lonlat = TRUE`) and `y` `NULL`;
  or `x` an sf or sfc object of POINT geometries, whose coordinate
  reference system must agree with `lonlat`. An unnamed two-column
  matrix is read as `x` then `y`; other columns are never taken by
  position. Vectors given separately (`x`, `y`, and `w` or `t` where a
  function takes them) are matched by name when every one is named and
  no name repeats, as elsewhere in the package; otherwise by position.
  No missing or infinite values.

- t:

  The times: numbers (days, say), or `Date` or `POSIXct` (taken in
  days), one per incident.

- ds:

  The upper edges of the distance bands: positive and increasing, in the
  unit of the coordinates (metres when `lonlat = TRUE`).

- dt:

  The upper edges of the time bands: positive and increasing, in the
  unit of `t`.

- nsim:

  The number of permutations: a whole number from 0 (none).

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

- lonlat:

  `FALSE` (the default) when `x` and `y` are planar coordinates in one
  length unit; `TRUE` when they are longitude and latitude in degrees,
  which are then projected to metres by
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  about their own centre, and centres are returned in degrees. It is
  never guessed.

## Value

A list of class `knxnr`: `tab`, a data frame with one row per cell,
distance bands varying fastest: `dlo`, `dhi`, `tlo`, `thi`, `obs`
(pairs), `exp`, `ratio` (`obs / exp`, `NaN` when nothing is expected),
`p` (`(1 + #{obs* >= obs}) / (nsim + 1)`, `NaN` without permutations)
and `ppois`; `n` (incidents), `pairs` (pairs within the largest
distance) and `nsim`.

## References

Knox, E. G. (1964). The detection of space-time interactions. *Journal
of the Royal Statistical Society, Series C*, 13(1), 25-30.
[doi:10.2307/2985220](https://doi.org/10.2307/2985220)

Townsley, M., Homel, R. and Chaseling, J. (2003). Infectious burglaries:
a test of the near repeat hypothesis. *British Journal of Criminology*,
43(3), 615-633.

## Examples

``` r
set.seed(1)
p <- cbind(x = runif(150, 0, 2000), y = runif(150, 0, 2000))
k <- knxnr(p, t = runif(150, 0, 365), ds = c(200, 400), dt = c(14, 28),
           nsim = 99, seed = 2)
k
#> knxnr: 150 incidents, 1322 close pairs, 99 permutations
#>  dlo dhi tlo thi obs      exp     ratio    p      ppois
#>    0 200   0  14  31 25.41369 1.2198149 0.12 0.15624997
#>  200 400   0  14  94 76.08752 1.2354195 0.02 0.02588842
#>    0 200  14  28  29 23.99195 1.2087390 0.16 0.17702659
#>  200 400  14  28  56 71.83087 0.7796091 0.99 0.97655462
```
