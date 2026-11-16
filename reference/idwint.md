# Inverse-distance-weighted interpolation, and zonal means

Carries the values of monitors to target points as the mean of the
nearest monitors weighted by their distance to the power `-p`, and
averages the targets of each zone. Exposure papers assign monitor
readings to places this way: Rahman et al. (2023, sec. 2.3) take
"inverse distance squared weighting of observations from up to four
nearby (\<50 km) monitoring locations" for census tract centroids
(`p = 2`, `k = 4`, `rad = 50`, `lonlat = TRUE`); Wongchanapai et al.
(2025) interpolate a surface of annual PM2.5 "using the inverse distance
weight method" and take "zonal statistics" per municipality (targets on
a grid over the state, `zone` the municipality of each grid point; they
used ArcGIS Pro and state neither the power nor the neighbours, so give
both); Wesselbaum (2022) matches each crime "to the measures of the
nearest weather, pollution station" (`k = 1`).

## Usage

``` r
idwint(
  mx,
  my,
  mv,
  tx,
  ty,
  p = 2,
  k = 0,
  rad = Inf,
  lonlat = FALSE,
  zone = NULL
)
```

## Arguments

- mx, my:

  The coordinates of the monitors.

- mv:

  The value of each monitor (`NA` allowed).

- tx, ty:

  The coordinates of the targets.

- p:

  The power of the distance, from 0.

- k:

  The number of nearest monitors used, a whole number from 0 (0: all
  within `rad`).

- rad:

  The radius: a positive number or `Inf`.

- lonlat:

  `TRUE` for longitude and latitude in degrees (distances in kilometres
  on the sphere); else planar.

- zone:

  `NULL` or the zone of each target (labels).

## Value

Without `zone`, a data frame, one row per target: `val` and `used`, the
monitors used. With `zone`, one row per zone (sorted): `zone`, `val`
(the zonal mean) and `n`, the targets with a value.

## Details

A target at a monitor's location takes that monitor's value (the mean of
the monitors there); else the `k` nearest monitors within `rad`
(`k = 0`: all within `rad`), ties in distance broken by the order of the
monitors, weigh `d^-p` each; `p = 0` gives their plain mean. Monitors
with a missing value are left out. A target with no monitor within `rad`
is `NA`. With `zone`, the zonal mean is the mean of the zone's targets
that have a value: with targets on a regular grid, the mean of the
interpolated surface over the zone.

## References

Rahman, M. M., Lorenzo, M., Ban-Weiss, G., Hasan, Z., Azzouz, M., Eckel,
S. P., Conti, D. V., Lurmann, F., Schlaerth, H., Johnston, J., Ko, J.,
Palinkas, L., Hurlburt, M., Silva, S., Gauderman, W. J., McConnell, R.
and Garcia, E. (2023). Ambient temperature and air pollution
associations with suicide and homicide mortality in California: a
statewide case-crossover study. *Science of the Total Environment*, 874,
162462.
[doi:10.1016/j.scitotenv.2023.162462](https://doi.org/10.1016/j.scitotenv.2023.162462)

Wongchanapai, P., Green, R., Ku, C. S., Peltier, M. R. and Ranjan, S.
(2025). Air pollution and intimate partner violence in New Jersey:
examining the evidence. *Psychology of Violence*, 15(6), 682-688.
[doi:10.1037/vio0000653](https://doi.org/10.1037/vio0000653)

Wesselbaum, D. (2022). Violent crimes and homicide in New York City: the
role of weather and pollution. *Journal of Forensic and Legal Medicine*,
91, 102430.
[doi:10.1016/j.jflm.2022.102430](https://doi.org/10.1016/j.jflm.2022.102430)

## See also

[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`ccross()`](https://rootcoder007.github.io/crimlab/reference/ccross.md)

## Examples

``` r
m <- data.frame(x = c(0, 10, 0), y = c(0, 0, 10), pm = c(8, 12, 10))
idwint(m$x, m$y, m$pm, c(5, 0), c(0, 0))
#>   val used
#> 1  10    3
#> 2   8    1
# the nearest monitor only
idwint(m$x, m$y, m$pm, c(9, 1), c(0, 1), k = 1)
#>   val used
#> 1  12    1
#> 2   8    1
# a grid over two zones, and the zonal means
g <- expand.grid(x = 0:10, y = 0:10)
idwint(m$x, m$y, m$pm, g$x, g$y, zone = ifelse(g$x < 5, "west", "east"))
#>   zone      val  n
#> 1 east 10.69539 66
#> 2 west  9.33759 55
```
