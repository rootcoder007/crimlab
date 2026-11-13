# Count regression: rate ratios between groups

Fits a log-linear model of a count on a group, with adjustment columns,
an exposure and clusters, by Poisson and negative binomial maximum
likelihood, and reports the incidence rate ratio (IRR) of each group
against the reference with its interval. The same model serves counts of
people, placements, days or deaths in the cells of a published table and
counts per person or per area.

## Usage

``` r
cntreg(
  y,
  grp = NULL,
  adj = NULL,
  expo = NULL,
  clus = NULL,
  ref = NULL,
  gee = FALSE,
  cadj = TRUE,
  level = 0.95
)
```

## Arguments

- y:

  The counts: finite numbers from 0, not all 0; or a `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md),
  whose elements fill the arguments by name (arguments given here
  override them).

- grp:

  The group of each count: labels (a factor by its levels, text or
  numbers by their sorted values), at least two groups.

- adj:

  `NULL` or adjustment columns: a data frame or list of columns, a
  numeric matrix, or one vector, one value per count. Factor, text and
  logical columns enter as indicators of each level but the first; give
  a number such as a year as a factor for one effect per value.

- expo:

  `NULL` or the exposure of each count (person-days, population):
  positive numbers; its logarithm is the offset.

- clus:

  `NULL` or the cluster of each count (institution, region, year):
  labels, at least two clusters.

- ref:

  The reference group: a label or a position; the first group by
  default.

- gee:

  `TRUE` to fit the exchangeable GEE as well (needs `clus`).

- cadj:

  `TRUE` for the `G / (G - 1)` factor on the cluster-robust variances.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `cntreg`, a list: `irr`, a data frame with one row
per fit and group (`fit`, `term`, `irr`, `lo`, `hi`, `p`, `sep`);
`coef`, every coefficient of every fit on the log scale (`fit`, `term`,
`est`, `se`, `z`, `p`); `fit`, one row per model (`fit`, `loglik`,
`aic`, `theta`, `setheta`, `alpha`, `iter`, `conv`, `rank`); `ref`, the
reference group; `n`, the number of counts; `ncl`, the number of
clusters (0 without `clus`); and `level`.

## Details

The Poisson fit takes the steps of
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) (iteratively
reweighted least squares from `y + 0.1`, stopping when the relative
change of the deviance is below 1e-8), and the negative binomial fit
(variance `mu + mu^2 / theta`) those of
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html): `theta` by
maximum likelihood and the coefficients in turn. A column of the design
that the columns before it explain (an adjustment level the group
determines, say) is aliased and dropped, as R does: its coefficients are
`NA`.

With `clus`, each fit also gets cluster-robust (sandwich) standard
errors, rows `"poiscl"` and `"nbcl"`:
`(X'WX)^-1 (sum_g U_g U_g') (X'WX)^-1` over the clusters' scores `U_g`,
times `G / (G - 1)` for `G` clusters when `cadj` (the default of
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html);
`cadj = FALSE` gives the independence GEE variance of 'geepack'). With
`gee = TRUE` the Poisson model is also fitted by generalised estimating
equations with an exchangeable working correlation (row `"gee"`), as
`geepack::geeglm(corstr = "exchangeable")` estimates it.

Intervals are Wald intervals on the log scale, `exp(est +- z se)`, with
`z` the normal quantile for `level`. `sep` flags a ratio whose group or
reference has no counts at all: the estimate then runs off towards 0 or
infinity and its interval means nothing.

Counts need not be whole numbers (a median duration, say): the
estimating equations are those of the Poisson pseudo-likelihood, and the
log-likelihood uses `lgamma(y + 1)`, where
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) reports an infinite
AIC.

## References

Venables, W. N. and Ripley, B. D. (2002). *Modern Applied Statistics
with S*, 4th ed. Springer.
[doi:10.1007/978-0-387-21706-2](https://doi.org/10.1007/978-0-387-21706-2)

Liang, K.-Y. and Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika*, 73(1), 13-22.
[doi:10.1093/biomet/73.1.13](https://doi.org/10.1093/biomet/73.1.13)

Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources*, 50(2), 317-372.
[doi:10.3368/jhr.50.2.317](https://doi.org/10.3368/jhr.50.2.317)

## Examples

``` r
# placements by group in eight institutions over three years
d <- data.frame(
  n = c(12, 30, 8, 25, 14, 41, 9, 22, 11, 35, 7, 28),
  female = rep(c("no", "yes"), 6),
  year = factor(rep(c(2021, 2022, 2023), each = 4)),
  inst = rep(c("A", "B", "C", "D"), each = 2, length.out = 12)
)
cntreg(d$n, d$female, adj = d["year"], clus = d$inst)
#> cntreg: 12 counts, 4 clusters; reference no
#> rate ratios (95% intervals):
#>     fit term      irr       lo       hi            p   sep
#>    pois  yes 2.967213 2.219882 3.966137 2.036053e-13 FALSE
#>      nb  yes 2.967220 2.211509 3.981170 4.108785e-13 FALSE
#>  poiscl  yes 2.967213 2.610266 3.372972 4.107025e-62 FALSE
#>    nbcl  yes 2.967220 2.610217 3.373050 4.296535e-62 FALSE
#> negative binomial theta 581.5672; AIC Poisson 74.68582, negative binomial 76.67518
# rates: a count per person-days of exposure
cntreg(c(5, 12), c("a", "b"), expo = c(1000, 1500))
#> Warning: cntreg: the negative binomial fit failed (theta has no finite estimate); its rows are NA.
#> cntreg: 2 counts; reference a
#> rate ratios (95% intervals):
#>   fit term irr       lo       hi         p   sep
#>  pois    b 1.6 0.563681 4.541576 0.3772454 FALSE
#>    nb    b  NA       NA       NA        NA FALSE
```
