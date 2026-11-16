# Smoke instruments from plume data

Builds the wildfire-smoke instruments of the recent papers from a daily
smoke-plume record by place: the presence of a plume; the presence
adjusted to the plumes that reached the ground; plume density classes;
smoke-attributable fine particles; and counts by month.

## Usage

``` r
smkiv(
  data,
  plume,
  unit,
  day,
  pm = NULL,
  adj = c("raw", "own", "nbr"),
  nbr = NULL,
  z = 1.64,
  seas = NULL,
  dens = NULL,
  cf = c("mean", "median", "childs"),
  by = NULL
)
```

## Arguments

- data:

  A data frame, one row per place and day.

- plume:

  The name of the plume column (0 or 1).

- unit, day:

  The names of the place and day columns.

- pm:

  `NULL` or the name of the PM2.5 column (needed for `"own"`, `"nbr"`
  and `spm`).

- adj:

  `"raw"`, `"own"` or `"nbr"`.

- nbr:

  For `"nbr"`, the name of the column naming each row's nearest
  neighbouring place (a value of `unit`).

- z:

  The number of standard deviations above the background.

- seas:

  `NULL` (meteorological seasons of a date `day`: dates, or ISO date
  text such as `"2024-06-01"`; a `day` of numbers has no seasons) or the
  name of a season column.

- dens:

  `NULL` or the name of a density-class column (`NA` or `""` for no
  plume).

- cf:

  `"mean"`, `"median"` or `"childs"`: the background of `spm`.

- by:

  `NULL` or the name of a grouping column (place-month) to aggregate to.

## Value

A data frame, one row per row of `data` (or per `by` group): `smoke` (0
or 1, or the count of smoke days per group), `bg` (the background, with
`pm`), `spm` (with `pm`), and one indicator (or count) per density
class, `dens:<class>`.

## Details

The choices and their sources:

- `adj = "raw"`: the plume indicator as given. Kircheis (2026, p. 4)
  instruments daily PM2.5 by "a daily binary indicator of predicted
  ground-level wildfire smoke presence" from a machine-learning product;
  Lee and Lee (2026, p. 6) use the same smoke days and, as a robustness
  check, NOAA's Hazard Mapping System plumes alone.

- `adj = "own"`: a plume day counts only when the place's PM2.5 is more
  than `z` standard deviations above its background, the mean of its
  days without a plume in the same season (Burkhardt et al. 2019,
  footnote 13: the own-county adjusted variable, 1.64 standard
  deviations, the 95th percentile under normality).

- `adj = "nbr"`: the same test made with the plume and the PM2.5 of the
  nearest neighbouring place (`nbr`), Burkhardt et al.'s preferred
  neighbour-adjusted variable (p. 11-12), which avoids using the place's
  own PM2.5 in its instrument.

- `seas`: the season of the background. Burkhardt et al. say "three
  month seasonal mean" without naming the months; the default is the
  meteorological seasons (December to February, ...) from a date column,
  pooled over years; give your own column to choose otherwise.

- `dens`: density classes (HMS "light", "medium", "heavy") become one
  indicator each. Lee and Lee (2026, footnote 4) warn that HMS density
  before 2018 has large measurement error, and use presence only.

- `spm`: smoke-attributable PM2.5, PM2.5 less the place-season
  background (`cf = "mean"`, Burkhardt et al.'s background, or
  `"median"`), on adjusted smoke days and at least 0; 0 on other days.
  This counterfactual is crimlab's construction from Burkhardt et al.'s
  background. `cf = "childs"` is instead the definition of Childs et al.
  (2022, equations 1-3, p. 13609): the background is the median PM2.5 of
  the same place and calendar month over the year and the years either
  side (three years) on days without smoke (`plume` 0), and smoke PM2.5
  is PM2.5 above it on smoke days (`plume` 1), at least 0, and 0 on
  other days; a place-month-year with no smoke-free days in its window
  gets `NA`, with a warning. It needs `day` as dates and ignores `seas`
  and `adj` (their smoke days are the plume column as given). Their
  machine-learned gridded smoke PM2.5 is a data product to be supplied
  as a column, not rebuilt here.

- `by`: with a group (place and month), the counts of smoke days and the
  means of the rest per group, Burkhardt et al.'s county-month
  instrument, "the number of days a county is affected by a fire smoke
  plume in a given month" (p. 10).

Bao and Da (2026, sec. 3.3) keep only "Medium and Heavy plumes" because
"including Light plumes badly overstates surface-smoke days" (give the
density classes in `dens` and use their indicators), and instead of a
plume overhead count the days a county lies downwind of a plume
([`dwnexp()`](https://rootcoder007.github.io/crimlab/reference/dwnexp.md)).

The instrument must move the exposure (a first stage, see
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md))
and change the outcome only through it: smoke may also keep people
indoors or darken the day (Lee and Lee 2026 discuss both), which the
exclusion restriction rules out by assumption.

## References

Childs, M. L., Li, J., Wen, J., Heft-Neal, S., Driscoll, A., Wang, S.,
Gould, C. F., Qiu, M., Burney, J. and Burke, M. (2022). Daily
local-level estimates of ambient wildfire smoke PM2.5 for the contiguous
US. *Environmental Science and Technology*, 56(19), 13607-13621.
[doi:10.1021/acs.est.2c02934](https://doi.org/10.1021/acs.est.2c02934)

Burkhardt, J., Bayham, J., Wilson, A., Berman, J. D., O'Dell, K., Ford,
B., Fischer, E. V. and Pierce, J. R. (2019). The relationship between
monthly air pollution and violent crime across the United States.
*Journal of Environmental Economics and Policy*, 9(2), 188-205.
[doi:10.1080/21606544.2019.1630014](https://doi.org/10.1080/21606544.2019.1630014)

Kircheis, L. (2026). Wildfire smoke increases assaults: evidence from
Seattle. *Environmental Research Letters*, 21, 044020.
[doi:10.1088/1748-9326/ae436c](https://doi.org/10.1088/1748-9326/ae436c)

Lee, G. and Lee, S. (2026). The impact of wildfire smoke exposure on
crime. *Environmental and Resource Economics*, 89, 5.
[doi:10.1007/s10640-025-01053-2](https://doi.org/10.1007/s10640-025-01053-2)

Bao, J. and Da, Y. (2026). Wildfire smoke and crime: evidence from
California. Working paper, Renmin University of China.

## See also

[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md),
[`exprun()`](https://rootcoder007.github.io/crimlab/reference/exprun.md),
[`dwnexp()`](https://rootcoder007.github.io/crimlab/reference/dwnexp.md)

## Examples

``` r
# simulated plumes and PM2.5 in two places over 120 days
set.seed(6)
d <- data.frame(place = rep(c("a", "b"), each = 120),
                day = rep(as.Date("2024-06-01") + 0:119, 2))
d$plume <- rbinom(240, 1, 0.2)
d$pm <- rnorm(240, 8, 2) + d$plume * rbinom(240, 1, 0.5) * 15
d$near <- ifelse(d$place == "a", "b", "a")
s <- smkiv(d, "plume", "place", "day", pm = "pm", adj = "own")
table(raw = d$plume, adjusted = s$smoke)
#>    adjusted
#> raw   0   1
#>   0 186   0
#>   1  25  29
head(smkiv(d, "plume", "place", "day", pm = "pm", adj = "nbr", nbr = "near"))
#>   smoke       bg spm
#> 1     0 7.978119   0
#> 2     0 7.978119   0
#> 3     0 7.978119   0
#> 4     0 7.978119   0
#> 5     0 7.978119   0
#> 6     0 7.978119   0
```
