# Downwind exposure: places downwind of sources of pollution

Whether each place lies downwind of sources of pollution in each period,
and how far off the wind's line it lies: the instruments of the papers
that let the wind carry smoke from fires that move from day to day (Bao
and Da 2026) or from a fixed set of sources (Ayesh 2023), for the
instrumental-variable estimators
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)
and
[`expgmm()`](https://rootcoder007.github.io/crimlab/reference/expgmm.md).

## Usage

``` r
dwnexp(
  sx,
  sy,
  sper,
  sdir,
  rx,
  ry,
  rid = NULL,
  rule = "sector",
  half = 30,
  rad = Inf,
  from = TRUE,
  lonlat = FALSE,
  wrap = TRUE,
  dwt = "mult"
)
```

## Arguments

- sx, sy:

  The coordinates of the sources.

- sper:

  `NULL` (one period) or the period of each source (labels).

- sdir:

  The wind direction at each source, degrees clockwise from north.

- rx, ry:

  The coordinates of the places: points, or the vertices of their
  polygons for `rule = "wedge"`.

- rid:

  `NULL` (each row a place) or the place of each row of `rx` (labels);
  for `rule = "wedge"`, the vertices of a place in order.

- rule:

  `"sector"`, `"wedge"` or `"angle"`.

- half:

  The half-angle of the sector, degrees from 0 to 180.

- rad:

  The radius of the sector: a positive number or `Inf`.

- from:

  `TRUE` when `sdir` is where the wind blows from.

- lonlat:

  `TRUE` for longitude and latitude in degrees (latitudes from -90 to
  90).

- wrap:

  `TRUE` for the angle about the circle (`rule = "angle"`).

- dwt:

  `"mult"`, `"div"` or `"none"`: the distance weight of the angle
  (`rule = "angle"`).

## Value

A data frame, one row per place and period (places fastest): `id`,
`per`, and `n` and `any` for `"sector"` and `"wedge"`, or `angle`,
`dist` and `z` for `"angle"` (`NA` where a period has no source at a
distance).

## Details

The sources are rows (`sx`, `sy`) with a period each (`sper`: a day, a
plume-day's date, a year) and the wind there (`sdir`, degrees clockwise
from north, where the wind blows from with `from = TRUE`: `sdir + 180`
is then the direction it blows to). The places are points (`rx`, `ry`)
or, for `rule = "wedge"`, polygons whose vertices are the rows of each
`rid` in order. Three rules:

- `"sector"`: a place is downwind of a source when the bearing from the
  source to the place is within `half` degrees of the direction the wind
  blows to, and the place within `rad` of the source (a place at the
  source counts). This is the point ("centroid") instrument of Bao and
  Da (2026, sec. 3.3), each county collapsed to one point and "the
  great-circle bearing from the plume anchor to that point" within 30
  degrees either side of the mean wind to-direction (`half = 30`,
  `rad = Inf`).

- `"wedge"`: a place's polygon meets the wedge of half-angle `half` and
  radius `rad` about that direction from the source: their baseline
  instrument, "a wedge-shaped transport corridor" with a half-angle of
  30 degrees and a radius of 1,000 km, a county downwind "if its census
  polygon intersects this wedge" (sec. 3.3). Planar coordinates only
  (projected, `rad` in their unit). The polygon meets the wedge when a
  vertex lies in it, the apex lies in the polygon, or an edge crosses a
  side or the arc.

- `"angle"`: per period, the angle between the bearing from the mean
  point of the sources to the place and the direction of their mean wind
  vector, and the mean distance from the place to the sources (those at
  its own location left out): Ayesh's (2023, sec. 4) instrument, "the
  absolute value of `lat_j - DirT_t`" with `DirT_t` the "average annual
  wind direction ... from the twelve rice growing districts", weighted
  by "the mean distance from the rice growing districts" (for a source
  district "the mean distance from eleven other rice growing
  districts", n. 23).

For `"sector"` and `"wedge"`, `n` counts the period's sources that reach
the place and `any` is 1 when one does; Bao and Da's instrument is the
sum of `any` over the days of a county-month (their eq 1), which
[`stats::aggregate()`](https://rdrr.io/r/stats/aggregate.html) gives
(see the examples). With several plumes a day, give each its own row
with that day's period.

The choices the papers leave open, as arguments:

- `from`: Bao and Da convert the "from-direction" to a "to-direction" by
  adding 180 degrees; Ayesh subtracts the wind direction from the
  azimuth without saying which way the direction is recorded.
  `from = TRUE` (the default) adds 180 degrees, so that the angle is 0
  for a place straight downwind; `from = FALSE` uses `sdir` as given.

- `wrap`: the angle about the circle, at most 180 (the default), or the
  plain absolute difference of the two bearings in `[0, 360)`, Ayesh's
  formula as written.

- `dwt`: Ayesh does not say how the mean distance weights the angle:
  `"mult"` (the default) gives `z = angle * dist`, `"div"`
  `angle / dist`, `"none"` the angle alone.

- `lonlat`: with `TRUE`, `x` and `y` are longitude and latitude in
  degrees, bearings are initial great-circle bearings and distances
  kilometres (haversine, radius 6371.0088 km); else planar, distances in
  the coordinates' unit.

The direction of the wind at a source is the user's: Bao and Da take the
speed-weighted circular mean of the to-direction over the wind grid
cells inside each plume
([`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md)
averages directions the same way), and keep only plumes of medium and
heavy density
([`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md)
reads those classes). As instruments these are valid only if the wind
carries the outcome nowhere but through the exposure.

## References

Bao, J. and Da, Y. (2026). Wildfire smoke and crime: evidence from
California. Working paper, Renmin University of China.

Ayesh, A. (2023). Burned agricultural biomass, air pollution and crime.
*Journal of Environmental Economics and Management*, 122, 102887.
[doi:10.1016/j.jeem.2023.102887](https://doi.org/10.1016/j.jeem.2023.102887)

## See also

[`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md),
[`smkiv()`](https://rootcoder007.github.io/crimlab/reference/smkiv.md),
[`updwn()`](https://rootcoder007.github.io/crimlab/reference/updwn.md),
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)

## Examples

``` r
# two plumes on day 1, one on day 2; four county points (planar, km)
s <- data.frame(x = c(0, 50, 0), y = c(0, 0, 0), day = c(1, 1, 2),
                from = c(180, 270, 0))
p <- data.frame(x = c(0, 100, 0, -80), y = c(100, 10, -60, 0))
r <- dwnexp(s$x, s$y, s$day, s$from, p$x, p$y)
r
#>   id per n any
#> 1  1   1 1   1
#> 2  2   1 1   1
#> 3  3   1 0   0
#> 4  4   1 0   0
#> 5  1   2 0   0
#> 6  2   2 0   0
#> 7  3   2 1   1
#> 8  4   2 0   0
# the county-month count of downwind days (Bao and Da 2026, eq 1)
r$month <- 1
aggregate(any ~ id + month, r, sum)
#>   id month any
#> 1  1     1   1
#> 2  2     1   1
#> 3  3     1   1
#> 4  4     1   0

# Ayesh's instrument: the angle off the mean wind of the source
# districts, times the mean distance to them
dwnexp(c(0, 10), c(0, 0), c(2020, 2020), c(240, 250), c(30, -30, 5),
       c(20, -20, 0), rule = "angle")
#>   id  per     angle     dist         z
#> 1  1 2020  13.65981 32.16989  439.4346
#> 2  2 2020 175.25512 40.38844 7078.2802
#> 3  3 2020  65.00000  5.00000  325.0000

# the wedge against a square county
sq <- data.frame(id = 1, x = c(-10, 10, 10, -10), y = c(40, 40, 60, 60))
dwnexp(0, 0, NULL, 180, sq$x, sq$y, sq$id, rule = "wedge", rad = 100)
#>   id per n any
#> 1  1   1 1   1
```
