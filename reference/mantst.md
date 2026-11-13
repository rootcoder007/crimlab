# Mantel's space-time test

The Pearson correlation, over all pairs of events, between their spatial
dissimilarity (the distance, or `1 / (distance + cs)` when `cs > 0`,
Mantel's reciprocal transform) and their temporal one (`|t_i - t_j|`, or
`1 / (|t_i - t_j| + ct)`), with a permutation test that permutes the
times over the fixed locations (Mantel 1967). Drawing distances and
intervals independently instead, as some descriptions do, does not keep
the dependence among the pairs of one distance matrix and gives no test
of stated size. Memory grows with the number of events, not pairs; time
with pairs times permutations.

## Usage

``` r
mantst(x, y = NULL, t, cs = 0, ct = 0, nsim = 999, seed = NULL, lonlat = FALSE)
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

  The time of each event: finite numbers (for dates,
  [`as.numeric()`](https://rdrr.io/r/base/numeric.html) of them, in days
  or seconds).

- cs, ct:

  The constants of the reciprocal transforms: numbers from 0 (0 uses the
  raw distance or interval).

- nsim:

  The number of permutations: one whole number from 0.

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

A named numeric vector: `r` (the correlation), `pairs`, and the
permutation p-values `pup` (`(1 + #{r_perm >= r}) / (nsim + 1)`, large
r: events close in space are close in time when both are raw distances),
`plo` (`r_perm <= r`) and `p2` (`|r_perm| >= |r|`); `NaN` without
permutations. 3 to 46340 events.

## References

Mantel, N. (1967). The detection of disease clustering and a generalized
regression approach. *Cancer Research*, 27(2), 209-220.

## Examples

``` r
set.seed(1)
x <- runif(60)
y <- runif(60)
mantst(x, y, t = x + rnorm(60, sd = 0.1), nsim = 199, seed = 1)
#>           r       pairs         pup         plo          p2 
#>    0.574889 1770.000000    0.005000    1.000000    0.005000 
```
