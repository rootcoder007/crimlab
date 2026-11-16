# Temperature inversions

The instrument of Bondy, Roth and Sager (2020, p. 566-567): "the
continuous difference between temperatures at the ground level and
temperatures at circa 750 m altitude (925 hectopascals (hPa) pressure
level)". On most days the air cools with height; in an inversion it is
warmer above than at the ground, and pollutants are trapped near the
ground.

## Usage

``` r
invsn(
  tlow,
  tup,
  sign = c("upper", "ground"),
  brk = NULL,
  by = NULL,
  rule = c("all", "any", "mean")
)
```

## Arguments

- tlow:

  Temperatures at the ground (numbers, `NA` allowed).

- tup:

  Temperatures aloft, at the same places and times, in the same unit.

- sign:

  `"upper"` or `"ground"`: which temperature comes first.

- brk:

  `NULL` or the breaks of the bins of the difference (as in
  [`cut()`](https://rdrr.io/r/base/cut.html), right-closed).

- by:

  `NULL` or the group of each reading (a day, say): one row per group,
  as below.

- rule:

  With `by`, when a group is an inversion: `"all"` (every reading with
  both temperatures inverted), `"any"` or `"mean"` (the mean difference
  above 0). Cautions from a review of 81 studies that use the instrument
  (Yamada 2026, section 7 and table 8, pp. 24-25), whose nine checks the
  help here follows. The instrument's strong first stages, mostly from
  China, do not carry over: in Sub-Saharan Africa and South Asia
  definitions drawn from published practice gave first stages that were
  strong yet second stages of opposite sign, and an inversion day moved
  PM2.5 by about one microgram per cubic metre, small beside the direct
  effects of the cold, still weather that comes with it, so the
  exclusion restriction fails through weather and season. The
  checks: (1) the reference layers must lie above ground everywhere in
  the study area; (2) the temperatures should be taken at local
  dawn; (3) report the share of zeros and the support of the instrument
  (the `diag` attribute of the result gives them); (4) estimate with
  both a pressure-level and a surface-layer definition and report
  both; (5) report the first stage in physical units relative to the
  mean exposure, not only its F (the `fstab` table of
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)); (6)
  where two definitions are available, test them against each other (an
  over-identification test across them; not computed here); (7) run a
  placebo outcome without a pollution pathway (any outcome in
  [`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)
  or
  [`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)); (8)
  check robustness to place-by-month fixed effects; (9) validate the
  exposure against ground monitors at the scale of the contrasts used.

## Value

A data frame: `diff`, the difference by `sign`; `inv`, 1 when `tup` is
above `tlow` and 0 otherwise (`NA` where a temperature is missing); and
with `brk`, `bin`, a factor. With `by`, one row per group, sorted: `by`,
`n` (the readings with both temperatures), `diff` (their mean
difference), `inv` (by `rule`; `NA` without a reading) and `bin`. The
attribute `"diag"` is a one-row data frame: `n`, the rows (groups) with
both temperatures; `zero`, the share of them without an inversion;
`min`, `max`, the support of `diff`.

## Details

The choices. The paper does not say which temperature is subtracted from
which: `sign = "upper"` (the default) gives the upper temperature less
the ground one, positive in an inversion; `sign = "ground"` the reverse,
as the paper's words run. `inv` is 1 in an inversion (upper warmer than
ground) whatever the sign; the indicator, rather than the continuous
difference, is the other common form of the instrument. With `brk`,
`bin` gives the strength in bins of the difference (by the chosen sign),
for an instrument that lets weak and strong inversions act differently.

From readings to days (`by`, `rule`). Batkeyev and DeRemer (2023, sec.
4.3, p. 297) compare a station at the city's elevation with one in its
high southeast corner, read every three hours, and define "a dummy which
equals 1 if within a day all low-elevation temperature readings are less
than high-elevation temperature readings" (`by` the day, `rule = "all"`,
the default), a demanding standard that still finds inversions on more
than 35% of their days; where evening readings were missing they use the
readings there are, as `"all"` does here. `"any"` and `"mean"` are
weaker rules, and `diff` (the mean difference) is the continuous form of
their appendix (sec. A3). Their other instrument, the share of a
period's readings with the wind from the mountains (from the southeast,
south-by-southeast or east-by-southeast, 101.25 to 168.75 degrees), set
to 0 for the districts far from the mountains, is
[`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md)
with `sect = c(101.25, 168.75)`, `by` the period and `share = TRUE`.

Cautions. Inversions cover large areas (the paper's satellite readings
are on a grid of about 110 by 70 km), so they vary little between places
in a city; Bondy, Roth and Sager use them at the borough and city level
and report a weaker first stage there than for wind (their table 5). As
an instrument the inversion must affect the outcome only through the
exposure; inversions come with cold, still, clear weather, which may
change what people do as well, so weather controls matter (see
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md)).

## References

Bondy, M., Roth, S. and Sager, L. (2020). Crime is in the air: the
contemporaneous relationship between air pollution and crime. *Journal
of the Association of Environmental and Resource Economists*, 7(3),
555-585. [doi:10.1086/707127](https://doi.org/10.1086/707127)

Yamada, H. (2026). When the air doesn't move the data: the thermal
inversion instrument in developing countries. Keio-IES Discussion Paper
DP2026-018, Keio University.

Batkeyev, B. and DeRemer, D. R. (2023). Mountains of evidence: the
effects of abnormal air pollution on crime. *Journal of Economic
Behavior and Organization*, 210, 288-319.
[doi:10.1016/j.jebo.2023.04.009](https://doi.org/10.1016/j.jebo.2023.04.009)

## See also

[`windiv()`](https://rootcoder007.github.io/crimlab/reference/windiv.md),
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)

## Examples

``` r
invsn(c(10, 5, 2), c(8, 6, 2))
#>   diff inv
#> 1   -2   0
#> 2    1   1
#> 3    0   0
# three-hourly readings, one row per day: inverted when every reading is
invsn(c(1, 2, 3, 1, 2, 3), c(4, 5, 6, 4, 1, 6), by = c(1, 1, 1, 2, 2, 2))
#>   by n     diff inv
#> 1  1 3 3.000000   1
#> 2  2 3 1.666667   0
invsn(c(10, 5, 2), c(8, 6, 2), sign = "ground", brk = c(-Inf, -1, 0, Inf))
#>   diff inv       bin
#> 1    2   0  (0, Inf]
#> 2   -1   1 (-Inf,-1]
#> 3    0   0    (-1,0]
# the share of days without an inversion and the support (check 3)
attr(invsn(c(10, 5, 2), c(8, 6, 2)), "diag")
#>   n      zero min max
#> 1 3 0.6666667  -2   1
```
