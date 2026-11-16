# Wooldridge's test for serial correlation in a panel

Whether the idiosyncratic errors of a linear panel model with unit
effects are serially correlated, by the test Wooldridge (2002) derives
and Drukker (2003) sets out and evaluates, the test Andreatta et al.
(2023, p. 1373) run on the residuals of their static models before
choosing a dynamic one. The first differences of `y` are regressed on
those of `x`, without a constant, which differencing removes (Drukker
2003, sec. 2); without serial correlation the residuals `e_it` of
consecutive differences correlate at -0.5, so `e_it` is regressed on
`e_i,t-1` with the variance clustered by unit, and the coefficient is
tested against -0.5: `F = ((coef + 0.5) / se)^2` on 1 and `G - 1`
degrees of freedom, `G` the units with a pair. A difference needs the
period before, so a gap removes two differences and the pairs that use
them. Drukker does not say whether the second regression has a constant;
`cons = FALSE` (the default) leaves it out, since the residuals of the
first have none to absorb.

## Usage

``` r
wsctst(y, x = NULL, unit, time, cons = FALSE, ssc = "cr1")
```

## Arguments

- y:

  The outcome of each row: finite numbers.

- x:

  `NULL` or the regressors: a numeric matrix or data frame, one row per
  row of `y`.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

- cons:

  `TRUE` for a constant in the regression of the residual on its lag.

- ssc:

  `"cr1"` or `"none"`: the small-sample factor of the variance.

## Value

An object of class `wsctst`, a list: `stat`, a named vector (`F`, `df1`,
`df2`, `p`, `coef` the coefficient of the lagged residual, `se`,
`pairs`, `clusters`, `ndiff` the differences); `coef`, the differenced
regression (`term`, `est`, `se`), `NULL` without `x`; and `ssc`.

## References

Drukker, D. M. (2003). Testing for serial correlation in linear
panel-data models. *The Stata Journal*, 3(2), 168-177.
[doi:10.1177/1536867X0300300206](https://doi.org/10.1177/1536867X0300300206)

Wooldridge, J. M. (2002). *Econometric Analysis of Cross Section and
Panel Data*. MIT Press.

Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
Digging into waste: an analysis of waste crime in the Italian provinces.
*Regional Studies*, 57(7), 1367-1379.
[doi:10.1080/00343404.2022.2135697](https://doi.org/10.1080/00343404.2022.2135697)

## Examples

``` r
set.seed(8)
d <- expand.grid(year = 1:8, prov = paste0("p", 1:40))
d$x <- rnorm(320)
d$e <- 0
for (i in which(d$year > 1)) d$e[i] <- 0.6 * d$e[i - 1] + rnorm(1)
d$y <- d$x + rnorm(40)[as.integer(substring(d$prov, 2))] + d$e
wsctst(d$y, d["x"], d$prov, d$year)
#> wsctst: Wooldridge test of serial correlation, F(1, 39) = 5.1486, p 0.02888
#> lagged residual -0.3814 (se 0.05228; -0.5 without correlation), 240 pairs
```
