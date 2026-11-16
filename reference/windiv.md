# Wind-direction instruments

Turns wind directions into indicators of the sector the wind blows from,
one sector left out, optionally one set per area so that the same wind
may act differently in different parts of a city; from hourly directions
it first makes one per day.

## Usage

``` r
windiv(
  dir,
  width = 90,
  start = 0,
  area = NULL,
  base = 1,
  sect = NULL,
  edges = c("halfopen", "paper"),
  by = NULL,
  avg = c("vector", "mode", "mean"),
  share = FALSE,
  spd = NULL
)
```

## Arguments

- dir:

  Wind directions in degrees (numbers; `NA` gives a row of `NA`, and is
  left out of a group's average).

- width:

  The width of a sector in degrees; it must divide 360 into at most
  67108864 (2^26) sectors.

- start:

  Where the first sector begins, in degrees.

- area:

  `NULL` or the area of each direction (labels; constant within a `by`
  group): one set of indicators per area, in sorted order of the labels.

- base:

  The sector left out, counted from 1.

- sect:

  `NULL`, or two numbers `c(from, to)` in degrees: then one indicator, 1
  when the wind blows from between `from` (included) and `to`
  (excluded), going clockwise; `width`, `start`, `base` and `edges` are
  not used.

- edges:

  `"halfopen"` or `"paper"`.

- by:

  `NULL`, or a group (day, or place and day) for each direction.

- avg:

  `"vector"`, `"mode"` or `"mean"`: how a group's directions become one.

- share:

  `TRUE` for the shares of each group's directions in each sector.

- spd:

  `NULL` or wind speeds, weights of the vector average.

## Value

A numeric matrix, one row per direction (or per `by` group, named by
it), columns named `"w<from>_<to>"`, with `":<area>"` added when `area`
is given. With `by`, the attribute `"dir"` holds the group's direction
(for `share`, its vector average).

## Details

The choices, and which paper made which:

Sectors (`width`, `start`, `edges`). Bondy, Roth and Sager (2020,
footnote 18) use quadrants "0-90, 91-180, 181-270" degrees with the
fourth left out, whole degrees closed at the top: `edges = "paper"`. The
default `"halfopen"` puts a boundary direction in the sector it opens,
`[lo, hi)`, which is the same for directions that are not whole numbers.
Bondy, Roth and Sager let each of five regions of London have its own
set (`area`), 15 instruments; their appendix also ranks four or eight
directions by each ward's long-run pollution (table A3), which is not
built here.

One sector (`sect`). Herrnstadt, Heyes, Muehlegger and Saberian (2019,
footnote 24) use one indicator, wind from between 225 and 285 degrees
(60 degrees); their table 5 widens it to 75 and 90 degrees without
saying where the wider sectors start, so give `sect` as you choose it
(for example `c(217.5, 292.5)` for 75 degrees about the same centre).

From hours to days (`by`, `avg`, `share`). With `by`, `dir` holds
several directions per day (hours) and each group gives one row.
`avg = "vector"` (the default) takes the direction of the mean wind
vector, `atan2(mean(sin), mean(cos))`, weighted by `spd` when given:
Herrnstadt and Muehlegger (2015, footnote 7) point out that 1 and 359
degrees average to north, not south. `avg = "mode"` takes the sector
most hours fall in (ties: the first), the "predominant wind direction"
of Bondy, Roth and Sager's figure A1. `avg = "mean"` is the plain mean
of the angles, wrong across north and given only for comparison.
`share = TRUE` instead gives the share of each group's hours in each
sector, Bondy, Roth and Sager's instrument ("the share of hours in the
24 hours of each day", footnote 18); `avg` is then not used.

Other papers that use the indicators, and their choices. Karamik and von
Graevenitz (2022, sec. 3.1) cut the circle into four 90-degree bins,
"North-East is equivalent to \[0 till 90)", with "the South-West ... the
most common wind direction in Germany" left out
(`width = 90, base = 3`), the direction averaged per region and day and
allowed to act differently in each region (`area`; their robustness
check narrows the bins to 60 degrees). Batkeyev and DeRemer (2023, sec.
4.3) use "the `[0,1]` share of weather observations within each time
period for which winds are blowing from the southeast,
south-by-southeast, or east-by-southeast", the three 22.5-degree points
from 101.25 to 168.75 degrees (`sect = c(101.25, 168.75)`, `by` the
period, `share = TRUE`), for the city's five southern districts and 0
for the northern ones. Singh and Visaria (2021, eq 4) assign each police
station the wind of the nearest grid point and use "eight wind direction
bins ... based on the octant from which the wind blows", interacted with
the nine administrative divisions of Bihar, with "wind from the east as
the base" (63 instruments); with octants centred on the compass points
that is `width = 45, start = 22.5, base = 2` and `area` the division
(the paper does not say where its octants start). Sarmiento (2020, eqs
7-8) enters "wind direction" in an hourly station-level first stage
without saying in what form; the sectors here are one form, and the
choice should be reported.

Directions are where the wind blows from, in degrees clockwise from
north; 360 is 0. The indicators are instruments only if the wind
direction affects the outcome through the exposure alone (see
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)).

## References

Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
contemporaneous relationship between air pollution and crime. *Journal
of the Association of Environmental and Resource Economists*, 7(3),
555-585. [doi:10.1086/707127](https://doi.org/10.1086/707127)

Herrnstadt, E., Heyes, A., Muehlegger, E. and Saberian, S. (2019). Air
pollution as a cause of violent crime. Working paper.

Herrnstadt, E. and Muehlegger, E. (2015). Air pollution and criminal
activity: evidence from Chicago microdata. NBER Working Paper 21787.
[doi:10.3386/w21787](https://doi.org/10.3386/w21787)

Karamik, Y. and von Graevenitz, K. (2022). Gone with the wind: the
effect of air pollution on crime. Evidence from Germany. ZEW Mannheim
working paper.
[doi:10.2139/ssrn.4130091](https://doi.org/10.2139/ssrn.4130091)

Batkeyev, B. and DeRemer, D. R. (2023). Mountains of evidence: the
effects of abnormal air pollution on crime. *Journal of Economic
Behavior and Organization*, 210, 288-319.
[doi:10.1016/j.jebo.2023.04.009](https://doi.org/10.1016/j.jebo.2023.04.009)

Sarmiento, L. (2020). I am innocent: hourly variations in air pollution
and crime behavior. DIW Berlin Discussion Paper.
[doi:10.2139/ssrn.3639084](https://doi.org/10.2139/ssrn.3639084)

Singh, T. P. and Visaria, S. (2021). Up in the air: air pollution and
crime. Evidence from India. Working paper, Hong Kong University of
Science and Technology.

## See also

[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md),
[`invsn()`](https://rootcoder007.github.io/crimlab/reference/invsn.md),
[`dwnexp()`](https://rootcoder007.github.io/crimlab/reference/dwnexp.md)

## Examples

``` r
windiv(c(10, 100, 200, 300))
#>      w90_180 w180_270 w270_0
#> [1,]       0        0      0
#> [2,]       1        0      0
#> [3,]       0        1      0
#> [4,]       0        0      1
# the paper's whole-degree quadrants: 90 is in the first
windiv(c(90, 91), edges = "paper")
#>      w90_180 w180_270 w270_0
#> [1,]       0        0      0
#> [2,]       1        0      0
# quadrants by region, as in London
windiv(c(10, 100, 200, 300), area = c("east", "east", "west", "west"))
#>      w90_180:east w180_270:east w270_0:east w90_180:west w180_270:west
#> [1,]            0             0           0            0             0
#> [2,]            1             0           0            0             0
#> [3,]            0             0           0            0             1
#> [4,]            0             0           0            0             0
#>      w270_0:west
#> [1,]           0
#> [2,]           0
#> [3,]           0
#> [4,]           1
# one westerly sector
windiv(c(230, 300), sect = c(225, 285))
#>      w225_285
#> [1,]        1
#> [2,]        0
# hours to days: 350 and 10 degrees average to north
windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2))
#>   w90_180 w180_270 w270_0
#> 1       0        0      0
#> 2       1        0      0
#> attr(,"dir")
#> [1] 360 110
windiv(c(350, 10, 100, 120), by = c(1, 1, 2, 2), share = TRUE)
#>   w90_180 w180_270 w270_0
#> 1       0        0    0.5
#> 2       1        0    0.0
#> attr(,"dir")
#> [1] 360 110
```
