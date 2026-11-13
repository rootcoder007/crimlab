# Trend in a series: rate ratio, Mann-Kendall and Theil-Sen

Measures the trend in a short series of counts (admissions, incidents or
placements per year) three ways. The Poisson log-linear model \\\log
\mu_t = b_0 + b_1 (x_t - \bar x) + \log e_t\\ gives the rate ratio per
unit of `x`, `exp(b1)` (0.85: a fall of 15% a period), with the exposure
`expo` as offset so that the trend is in the rate; when the Pearson
dispersion (over `n - 2`) is above 1.5 the interval and p-value are the
quasi-Poisson ones (standard error times the square root of the
dispersion, the t distribution on `n - 2` degrees of freedom). The fit
takes the steps of [`stats::glm()`](https://rdrr.io/r/stats/glm.html).
Mann-Kendall's rank test and the Theil-Sen slope (the median of the
pairwise slopes, with Sen's distribution-free interval) make no
distributional assumption; the p-value is exact over every ordering of
the observed values up to 8 periods (ties kept) and normal with the tie
and continuity corrections beyond. Negating the series negates S, so
rises and falls are treated alike (`proof("trndrr")`). `yoy` holds the
change on the period before, and with `seas` the totals by season are
tested for equality (chi-square).

## Usage

``` r
trndrr(y, x = NULL, expo = NULL, seas = NULL, alt = "both", level = 0.95)
```

## Arguments

- y:

  The counts, in period order: numbers from 0, not all 0, at least 3.

- x:

  `NULL` for the periods 1, 2, ..., or the periods: strictly increasing
  numbers (years, say).

- expo:

  `NULL` or the exposure of each period (population, person-days):
  positive numbers.

- seas:

  `NULL` or the season of each count (month, day of week): labels, two
  or more.

- alt:

  The alternative of the Mann-Kendall test: `"both"`, `"up"` or
  `"down"`.

- level:

  The confidence level: one number strictly between 0 and 1.

## Value

An object of class `trndrr`, a list: `rr`, the rate ratio row (`rr`,
`lo`, `hi`, `p`, the log slope `b` and its `se`, the dispersion `disp`,
`over` when the quasi-Poisson interval was used); `mk` (`s`, its
variance `var`, `tau`, `p`, `exact`); `sen` (`slope`, `int`, `lo`, `hi`,
`clamp`); `yoy` (`x`, `y`, the change `chg` and per cent change `pct` on
the period before); `seas`, `NULL` or a list of `tab` (`season`, `n`,
`expd`), `stat`, `df` and `p`; `fitted`, the Poisson fitted counts; `n`,
`alt` and `level`.

## Details

When every count falls in the first or the last period the likelihood
rises without bound: the rate ratio is 0 or `Inf`, with `NA` for the
rest of its row, and a warning. Sen's interval needs enough periods:
`clamp` is `TRUE` when an end had to be the extreme slope.

## References

Sen, P. K. (1968). Estimates of the regression coefficient based on
Kendall's tau. *Journal of the American Statistical Association*,
63(324), 1379-1389.
[doi:10.1080/01621459.1968.10480934](https://doi.org/10.1080/01621459.1968.10480934)

Hamed, K. H. and Rao, A. R. (1998). A modified Mann-Kendall trend test
for autocorrelated data. *Journal of Hydrology*, 204(1-4), 182-196.
[doi:10.1016/S0022-1694(97)00125-X](https://doi.org/10.1016/S0022-1694%2897%2900125-X)

## See also

[`stpcp()`](https://rootcoder007.github.io/crimlab/reference/stpcp.md)
for a step change,
[`xcser()`](https://rootcoder007.github.io/crimlab/reference/xcser.md)
to correlate two series.

## Examples

``` r
y <- c(402, 377, 390, 368, 331, 310, 305)
trndrr(y, 2017:2023)
#> trndrr: 7 periods
#> rate ratio per period 0.9523448 (95% 0.9337437 to 0.9713165), p 1.223921e-06
#> Mann-Kendall S -19, p 0.002777778; Theil-Sen slope -17.75
# a rate: counts over a growing population
trndrr(c(20, 25, 30, 33), expo = c(1000, 1500, 2500, 2600))$rr
#>          rr        lo       hi          p          b         se      disp  over
#> 1 0.8494552 0.7109716 1.014913 0.07234571 -0.1631601 0.09079901 0.4268193 FALSE
# equal totals by quarter?
trndrr(c(5, 9, 7, 4, 6, 10, 8, 5), seas = rep(1:4, 2))$seas
#> $tab
#>   season  n expd
#> 1      1 11 13.5
#> 2      2 19 13.5
#> 3      3 15 13.5
#> 4      4  9 13.5
#> 
#> $stat
#> [1] 4.37037
#> 
#> $df
#> [1] 3
#> 
#> $p
#> [1] 0.2241485
#> 
```
