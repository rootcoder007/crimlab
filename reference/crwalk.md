# Correlated random walk of a series of events

Treats the successive locations of a series (incidents attributed to one
offender, in time order) as a correlated random walk (Kareiva and
Shigesada 1983): steps of independent lengths and changes of heading
(turning angles) independent of each other and of the steps. The fit
gives the mean and mean squared step and the mean cosine and sine of the
turns; the expected squared displacement after `k` moves is Kareiva and
Shigesada's for symmetric turning (mean sine 0),
`k E l^2 + 2 (E l)^2 c / (1 - c) (k - (1 - c^k) / (1 - c))`, next to the
observed mean over every start. The walk is then continued `nsim` times
for `h` steps from the last location and heading, each turn and step
drawn independently and uniformly from the observed ones (an empirical
correlated random walk), and the simulated positions summarise where the
next events may fall. It describes movement; it does not identify
anyone, and a series with few events gives a rough walk.

## Usage

``` r
crwalk(x, y = NULL, h = 1, nsim = 999, seed = NULL, lonlat = FALSE)
```

## Arguments

- x, y:

  The points: `x` and `y` numeric vectors of equal length; or `x` a
  matrix or data frame with columns named `x` and `y` (or `X` and `Y`,
  as `sf::st_coordinates()` names them, or `lon` and `lat` with
  `lonlat = TRUE`) and `y` `NULL`; or `x` an sf or sfc object of POINT
  geometries, whose coordinate reference system must agree with
  `lonlat`. An unnamed two-column matrix is read as `x` then `y`; other
  columns are never taken by position. Vectors given separately (`x`,
  `y`, and `w` or `t` where a function takes them) are matched by name
  when every one is named and no name repeats, as elsewhere in the
  package; otherwise by position. No missing or infinite values.

- h:

  The number of steps ahead: a whole number from 1.

- nsim:

  The number of simulated continuations: a whole number from 1.

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

A list. `fit`, a named vector: `steps`, `ml` (mean step), `ml2` (mean
squared step), `c` and `s` (mean cosine and sine of the turning angles)
and `head` (the last heading, radians anticlockwise from east). `msd`, a
data frame with one row per step ahead `k`: `obs` (the observed mean
squared displacement over `k` moves, `NA` beyond the series), `exp` (its
expectation; it assumes `s` near 0), and `x`, `y`, `sdx`, `sdy` (the
mean and standard deviations of the simulated positions; the means in
degrees when `lonlat = TRUE`, the rest in metres). `sim`, the simulated
positions: matrices `x` and `y`, one row per walk and one column per
step (planar, metres when `lonlat = TRUE`).

## References

Kareiva, P. M. and Shigesada, N. (1983). Analyzing insect movement as a
correlated random walk. *Oecologia*, 56(2-3), 234-238.
[doi:10.1007/BF00379695](https://doi.org/10.1007/BF00379695)

## Examples

``` r
set.seed(1)
a <- cumsum(rnorm(20, 0, 0.5))
p <- cbind(x = cumsum(cos(a)), y = cumsum(sin(a)))
r <- crwalk(p, h = 3, nsim = 200, seed = 1)
r$fit
#>      steps         ml        ml2          c          s       head 
#> 19.0000000  1.0000000  1.0000000  0.8935287  0.1187582  1.9052388 
r$msd
#>   k      obs      exp        x        y       sdx       sdy
#> 1 1 1.000000 1.000000 12.12662 10.09144 0.3882481 0.1794948
#> 2 2 3.787057 3.787057 11.60345 10.74197 0.7778394 0.4737584
#> 3 3 8.165410 8.170902 11.07326 11.25329 1.1794886 0.8584427
```
