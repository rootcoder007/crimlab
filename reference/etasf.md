# ETAS: a self-exciting process with marks and the Omori law

Fits the temporal epidemic-type aftershock sequence (ETAS) model of
Ogata (1988) to event times with marks (a magnitude, the number of
victims, the severity of an incident): \$\$\lambda(t) = \mu + K
\sum\_{t_j \< t} e^{\alpha (m_j - m_0)} (t - t_j + c)^{-p},\$\$ so
bigger events trigger more, and the triggering decays as a power of the
time since (the modified Omori law). Without marks, or with constant
ones, \\\alpha\\ is fixed at 0 and the model is a Hawkes process with a
Lomax kernel (that of
[`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md)
with shape \\p - 1\\ and scale \\c\\). The log-likelihood over
`[0, hor]` is maximised with its analytic gradient by BFGS over the
logarithms of \\\mu\\, \\K\\, \\c\\ and \\p - 1\\ (and \\\alpha\\), from
two starts or `start`; standard errors from the numerical Hessian of the
gradient. `branch` is the branching ratio \\K \bar{e^{\alpha (m - m_0)}}
c^{1 - p} / (p - 1)\\, the mean number of events each event triggers. On
data with little power-law tail the likelihood rises towards the
exponential limit of the kernel (`c` and `p` growing together): the fit
then stops with `conv` `FALSE`, and
[`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md)'s
exponential kernel is the better model.

## Usage

``` r
etasf(times, hor, mag = NULL, m0 = NULL, start = NULL, maxit = 500)
```

## Arguments

- times:

  Event times: numbers sorted in `[0, hor]`, at least 10.

- hor:

  The end of the window: one positive number.

- mag:

  `NULL` or the mark of each event: finite numbers.

- m0:

  `NULL` (the least mark) or the reference mark: one number.

- start:

  `NULL` or the starting `mu`, `k`, `alpha`, `c` and `p`.

- maxit:

  The most BFGS steps per start: one whole number from 1.

## Value

An object of class `etasf`, a list: `par`, `se` and `cov` of (`mu`, `k`,
`alpha`, `c`, `p`; `alpha`'s are `NA` when it is fixed); `loglik`,
`aic`; `branch`; `conv`, `iter`; `afree`, whether `alpha` was estimated;
`n` and `hor`.

## References

Ogata, Y. (1988). Statistical models for earthquake occurrences and
residual analysis for point processes. *Journal of the American
Statistical Association*, 83(401), 9-27.
[doi:10.1080/01621459.1988.10478560](https://doi.org/10.1080/01621459.1988.10478560)

## See also

[`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md),
[`mvhwk()`](https://rootcoder007.github.io/crimlab/reference/mvhwk.md),
[`sthwk()`](https://rootcoder007.github.io/crimlab/reference/sthwk.md).

## Examples

``` r
x <- hwksim(c(log(0.3), 0.5, 0.6, 1), 800, kernel = "lomax", seed = 5)
etasf(x$t, 800)
#> etasf: 530 events in [0, 800]; no marks, alpha fixed at 0
#>             est         se
#> mu    0.3327782 0.08931391
#> k     0.4373506 0.87721720
#> alpha 0.0000000         NA
#> c     1.3824477 1.65692065
#> p     1.6791760 0.82290541
#> branching ratio 0.5167995; log-likelihood -724.0694, AIC 1456.139
```
