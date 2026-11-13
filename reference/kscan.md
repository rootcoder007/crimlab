# Space-time scan statistic of incidents

Kulldorff's space-time permutation scan (Kulldorff, Heffernan, Hartman,
Assuncao and Mostashari 2005), which needs the incidents alone, no
population: cylinders are circles of radius `r` about each distinct
incident location times windows `[s, s + tw)` that start at each
distinct incident time. A cylinder holding `c` incidents, of `n_C` in
its circle (at any time) and `n_W` in its window (anywhere), expects
`e = n_C n_W / n` when place and time are independent, and scores the
Poisson log-likelihood ratio
`c log(c / e) + (n - c) log((n - c) / (n - e))` when `c > e`, 0
otherwise (Kulldorff 1997). The most likely cluster has the largest
ratio; secondary clusters follow in order of their ratio, each the best
of its circle among those whose circles overlap none chosen before. Each
p-value compares the cluster's ratio with the largest ratio of the whole
scan with the times shuffled among the locations `nsim` times from
`seed`, so a secondary cluster's p-value is conservative. That the ratio
is never negative is proved in Lean: see `proof("kscan")`. The work
grows as the number of locations times `length(r)` times the distinct
times times `length(tw)`, for the scan and for each permutation.

## Usage

``` r
kscan(
  x,
  y = NULL,
  t,
  r,
  tw,
  minc = 2,
  top = 1,
  nsim = 999,
  seed = NULL,
  lonlat = FALSE
)
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

- r:

  The radii of the circles: positive and increasing, in the unit of the
  coordinates (metres when `lonlat = TRUE`).

- tw:

  The lengths of the windows: positive, in the unit of `t`.

- minc:

  The fewest incidents a cluster holds: a whole number from 1.

- top:

  The most clusters reported: a whole number from 1.

- nsim:

  The number of permutations: a whole number from 0.

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

A list of class `kscan`: `cl`, a data frame with one row per cluster,
most likely first: its centre `x`, `y` (in degrees when
`lonlat = TRUE`), `r`, `start` and `end` (the window `[start, end)`),
`n` (incidents in the cylinder), `exp`, `ncirc` and `nwin` (incidents in
its circle and in its window), `rr` (the relative risk inside against
outside), `llr` and `p` (`(1 + #{max* >= llr}) / (nsim + 1)`, `NaN`
without permutations); `memb`, each incident's cluster or `NA`; `n` and
`nsim`.

## References

Kulldorff, M. (1997). A spatial scan statistic. *Communications in
Statistics - Theory and Methods*, 26(6), 1481-1496.
[doi:10.1080/03610929708831995](https://doi.org/10.1080/03610929708831995)

Kulldorff, M., Heffernan, R., Hartman, J., Assuncao, R. and Mostashari,
F. (2005). A space-time permutation scan statistic for disease outbreak
detection. *PLoS Medicine*, 2(3), e59.
[doi:10.1371/journal.pmed.0020059](https://doi.org/10.1371/journal.pmed.0020059)

## Examples

``` r
set.seed(1)
p <- cbind(x = c(runif(150, 0, 20), rnorm(20, 5, 0.3)),
           y = c(runif(150, 0, 20), rnorm(20, 15, 0.3)))
t <- c(sample(0:99, 150, TRUE), sample(40:44, 20, TRUE))
s <- kscan(p, t = t, r = c(1, 3), tw = c(5, 20), top = 2, nsim = 99,
           seed = 1)
s
#> kscan: 2 clusters of 170 incidents, 99 permutations
#>         x        y r start end  n       exp ncirc nwin        rr       llr    p
#>  4.959446 15.23629 1    40  45 20 3.4117647    20   29  6.510345 19.636580 0.01
#>  9.544601 11.09801 3     4   9  4 0.3235294    11    5 12.637459  6.422689 0.12
```
