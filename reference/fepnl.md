# Two-way fixed effects with an interaction, for counts of offences

The regression Duan et al. (2026, eqs 6-8, pp. 9-11) use to check the
drivers of recorded pollution crimes: the count on continuous drivers
and the product of two of them, with unit and period fixed effects,
errors clustered by unit, standardised coefficients, variance inflation
factors, the simple slopes of one driver at values of the other (p. 23),
and the Jarque-Bera test of the residuals (p. 25).

## Usage

``` r
fepnl(
  y,
  x,
  unit,
  time,
  int = NULL,
  clus = NULL,
  at = NULL,
  ssc = "cr1",
  level = 0.95
)

# S3 method for class 'fepnl'
plot(x, ..., main = NULL, col = NULL, pal = NULL)
```

## Arguments

- y:

  The outcome of each row: finite numbers.

- x:

  The drivers: a numeric matrix or data frame, one row per row of `y`
  (or one vector); column names name the terms. For
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), a `fepnl`
  result.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers; each unit at most once per period.

- int:

  `NULL` or the names of two columns of `x` whose product is added: the
  slopes are those of the first at values of the second.

- clus:

  `NULL` (cluster by unit) or the cluster of each row (labels).

- at:

  `NULL` or the values of `int[2]` for the simple slopes.

- ssc:

  `"cr1"` or `"none"`: the small-sample factor of the variance.

- level:

  The confidence level of the intervals.

- ...:

  Passed to the base-graphics call.

- main:

  `NULL` (a default title) or one string.

- col:

  `NULL` (the palette) or colours.

- pal:

  `NULL` (the session palette, see
  [`plotopt()`](https://rootcoder007.github.io/crimlab/reference/blpal.md))
  or a palette name for
  [`blpal()`](https://rootcoder007.github.io/crimlab/reference/blpal.md).

## Value

An object of class `fepnl`, a list: `coef`, one row per term (`term`,
`est`, `se`, `lo`, `hi`, `t`, `p`, `beta`, `vif`); `slopes` (`at`,
`slope`, `se`, `lo`, `hi`, `t`, `p`), `NULL` without `int`; `fit`, one
row: `n`, `units`, `periods`, `ncl` (`G`), `k` (`K`), `r2`, `r2w`,
`rss`; `jb`, a named vector `skew`, `kurt`, `jb`, `p`; `vcov`; `resid`,
the residuals in the order of the rows; `ssc`; and `level`.

## Details

The model is `y_it = sum_j b_j x_itj + a_i` `+ m_t + e_it`, fitted by
least squares after the within transformation (the kernel of
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md),
so an unbalanced panel is fine); with `int = c("a", "b")` the column
`a:b`, the product of `a` and `b`, is added to `x`. The variance is CR1
by `clus` (`G / (G - 1) (n - 1) / (n - K)` times the sandwich, as in
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md));
`ssc = "none"` drops the factor. Duan et al. fitted the model with the
'linearmodels' PanelOLS estimator clustered by province and do not state
its small-sample factor. Tests and intervals use Student's t on `G - 1`
degrees of freedom.

`beta` is the standardised coefficient `b SD_x / SD_y` of Duan et al.'s
step 4 (p. 11), with the standard deviations of the columns as given
(the interaction column for its term). `vif` is the variance inflation
factor `1 / (1 - R2_j)` of each column of `x` (with the interaction)
regressed on a constant and the others over the rows used; Duan et al.
drop drivers with a VIF above 10 before the fit (p. 10) and do not say
whether their VIFs were pooled, as here, or within. The simple slope of
`int[1]` at a value `m` of `int[2]` is `b_a + b_ab m` with variance
`V_aa + m^2 V_ab,ab + 2 m V_a,ab`; by default at the mean and the mean
plus and minus one standard deviation of `int[2]`, Duan et al.'s points.
`r2` is `1 - RSS / TSS` of the full model (unit and period effects
counted as fitted), `r2w` the within R-squared. The Jarque-Bera
statistic is `n / 6 (S^2 + (K - 3)^2 / 4)` of the residuals, with
moments divided by `n`, on two degrees of freedom.

Duan et al.'s model is linear in a count. For a Poisson or negative
binomial model with the same effects see
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
(binary treatment) and
[`cntreg()`](https://rootcoder007.github.io/crimlab/reference/cntreg.md);
Andreatta et al. (2023) fit a dynamic panel of the logged rate by system
GMM, as
[`gmmpnl()`](https://rootcoder007.github.io/crimlab/reference/gmmpnl.md)
does. Melhim (2026) fits the annual crime rate of ten US states
(2000-2022) on PM2.5, poverty and unemployment by pooled least squares
and then with state and year fixed effects and errors clustered by
state, the model here without `int`; the estimate changes from
insignificant to positive once the state effects enter.

The plot draws the coefficients with their intervals and, with an
interaction, the simple slopes at each value of the moderator; it
returns what it drew.

## References

Duan, H., Nie, C., He, B., Chen, X., Duan, Z., Chen, S., Song, J. and
Gao, Z. (2026). A multi-method framework for decoding the spatiotemporal
patterns and driving factors of environmental pollution crimes in China.
*Humanities and Social Sciences Communications*.
[doi:10.1057/s41599-026-08766-9](https://doi.org/10.1057/s41599-026-08766-9)

Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
Digging into waste: an analysis of waste crime in the Italian provinces.
*Regional Studies*, 57(7), 1367-1379.
[doi:10.1080/00343404.2022.2135697](https://doi.org/10.1080/00343404.2022.2135697)

Melhim, S. (2026). From air quality to public safety: rethinking crime
dynamics and social costs of pollution in the United States. *Journal of
Environmental Health*, 89(3), 18-25.
[doi:10.70387/001c.171786](https://doi.org/10.70387/001c.171786)

## Examples

``` r
# 30 provinces over 8 years: cases on the case rate, the unemployment
# rate and their product, simulated
set.seed(2)
d <- expand.grid(year = 2016:2023, prov = paste0("p", 1:30))
d$cper <- runif(240, 0.005, 0.04)
d$ur <- runif(240, 0.01, 0.03)
d$eipa <- rnorm(240, 5)
d$cases <- 50 + 9000 * d$cper + 1500 * d$ur - 4e5 * d$cper * d$ur +
  rnorm(240, sd = 5)
r <- fepnl(d$cases, d[c("cper", "eipa", "ur")], d$prov, d$year,
           int = c("cper", "ur"))
r
#> fepnl: 240 rows, 30 units, 8 periods, 30 clusters, R2 0.9929, within R2 0.9916
#>     term           est           se            lo            hi           t
#>     cper  9.026466e+03  103.3399520  8.815113e+03  9.237820e+03  87.3473066
#>     eipa  1.888576e-01    0.3811923 -5.907681e-01  9.684834e-01   0.4954393
#>       ur  1.534870e+03  106.8914777  1.316252e+03  1.753488e+03  14.3591418
#>  cper:ur -3.997308e+05 4706.6920055 -4.093571e+05 -3.901046e+05 -84.9281898
#>             p         beta       vif
#>  1.128269e-36  1.785147290 11.771218
#>  6.240236e-01  0.003512598  1.001826
#>  1.030284e-14  0.174856822  5.533638
#>  2.539870e-36 -1.866182026 15.657058
#> simple slopes:
#>          at     slope       se        lo        hi         t            p
#>  0.01356247  3605.130 43.35037  3516.468  3693.791  83.16260 4.659674e-36
#>  0.01963659  1177.117 24.09165  1127.844  1226.390  48.85996 2.082308e-29
#>  0.02571071 -1250.895 30.26917 -1312.803 -1188.988 -41.32572 2.504726e-27
#> Jarque-Bera of the residuals 1.386, p 0.5002
r$slopes
#>           at     slope       se        lo        hi         t            p
#> 1 0.01356247  3605.130 43.35037  3516.468  3693.791  83.16260 4.659674e-36
#> 2 0.01963659  1177.117 24.09165  1127.844  1226.390  48.85996 2.082308e-29
#> 3 0.02571071 -1250.895 30.26917 -1312.803 -1188.988 -41.32572 2.504726e-27
```
