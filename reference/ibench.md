# Internal benchmarking of officers

Asks whether an officer stops more people of a group (or searches,
arrests, uses force more often) than other officers do in the same
circumstances. Officers work different places and hours, so each
officer's stops are compared not with all other stops but with the other
officers' stops weighted to look like this officer's: the propensity
`e(x)` that a stop with circumstances `x` (place, hour, day, reason) is
this officer's comes from a logistic regression on `x` (over all stops),
and each other officer's stop gets the weight \`e / (1

- e)`, the weights of the effect on the treated (Ridgeway and MacDonald 2009). The benchmark is the weighted share of the outcome among those stops, and the z statistic compares the officer's share with it, with the variance of a proportion for the officer and, for the benchmark, that of a proportion over its effective number of stops, `(sum
  w)^2 / sum
  w^2`. With many officers some will stand out by chance: `q`is the false discovery rate of Benjamini and Hochberg (1995) over the officers tested, and`flag`marks those with`q`below`fdr\`.

## Usage

``` r
ibench(y, off, x = NULL, minn = 20, fdr = 0.05, dr = FALSE, local = FALSE)
```

## Arguments

- y:

  The outcome of each stop: 0 or 1 (or `FALSE` or `TRUE`), such as
  whether the person stopped belongs to the group of interest.

- off:

  The officer of each stop: labels, two officers at least.

- x:

  `NULL` or the circumstances of each stop: a numeric matrix, a data
  frame or list of numeric columns, or one numeric vector, one row per
  stop (code a place or an hour as indicator columns).

- minn:

  The fewest stops of an officer to test: one whole number from 1.

- fdr:

  The false discovery rate at which an officer is flagged: one number
  strictly between 0 and 1.

- dr:

  `TRUE` for the doubly robust benchmark, `FALSE` for the propensity
  weights alone.

- local:

  `TRUE` to estimate the local false discovery rate and flag by it.

## Value

An object of class `ibench`, a list: `tab`, one row per officer (`off`,
`n` stops, `rate`, the officer's share of the outcome, `bench`, the
benchmark, `ess`, its effective number of stops, `z`, `p`, `q`, with
`local` `lfdr`, and `flag`); `fdr`; `minn`; `nx`, the number of columns
of `x`; `dr`; and `null`, `NULL` without `local`, otherwise the
empirical null's `delta`, `sigma` and `p0`.

## Details

Two options follow Ridgeway and MacDonald (2009) further:

- `dr = TRUE`, the doubly robust benchmark: the propensity weights are
  kept, and a logistic regression of the outcome on an indicator of the
  officer's stops and on `x`, weighted 1 for the officer's stops and
  `e / (1 - e)` for the others, adjusts for what the weights leave
  unbalanced. The benchmark is the mean of its fitted chances for the
  officer's stops with the indicator set to 0 (the same stops made by
  the others), and `z` is the officer's rate less the benchmark over the
  sandwich (robust) standard error of that difference from the weighted
  fit. The estimate is right when either the propensity model or the
  outcome model is. This is the propensity-weighted regression of the
  first author's fastDR package (which fits the propensity by boosting,
  where a logistic regression is used here).

- `local = TRUE`, the local false discovery rate of Efron (2004) for
  each tested officer: the chance that an officer with that `z` is one
  of the unremarkable majority. The density of all the officers' `z` is
  fitted by a Poisson regression on a natural spline of 7 df over a
  histogram of 120 breaks, the null density is a normal
  `N(delta, sigma^2)` fitted by maximum likelihood to the middle of the
  `z` values (the empirical null, which absorbs spread that the
  benchmark misses), and `lfdr = p0 f0(z) / f(z)`, as in the defaults of
  the package locfdr; `null` gives `delta`, `sigma` and `p0`. `flag`
  then marks the officers with `lfdr` below `fdr` (Efron suggests 0.2).
  The method is meant for many officers, hundreds rather than tens; when
  the null cannot be fitted the function refuses.

The comparison is only as good as `x`: circumstances that differ between
officers and are not in `x` are left in the difference. An officer with
fewer than `minn` stops, or whose propensity model does not converge
(circumstances that separate this officer's stops from all others), is
not tested (`NA`).

## References

Ridgeway, G. and MacDonald, J. M. (2009). Doubly robust internal
benchmarking and false discovery rates for detecting racial bias in
police stops. *Journal of the American Statistical Association*,
104(486), 661-668.
[doi:10.1198/jasa.2009.0034](https://doi.org/10.1198/jasa.2009.0034)

Benjamini, Y. and Hochberg, Y. (1995). Controlling the false discovery
rate: a practical and powerful approach to multiple testing. *Journal of
the Royal Statistical Society, Series B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Efron, B. (2004). Large-scale simultaneous hypothesis testing: the
choice of a null hypothesis. *Journal of the American Statistical
Association*, 99(465), 96-104.
[doi:10.1198/016214504000000089](https://doi.org/10.1198/016214504000000089)

## Examples

``` r
set.seed(2)
n <- 3000
off <- sample(sprintf("o%02d", 1:12), n, TRUE)
night <- rbinom(n, 1, ifelse(off %in% c("o01", "o02"), 0.7, 0.3))
y <- rbinom(n, 1, plogis(-1 + night + 0.8 * (off == "o03")))
ibench(y, off, data.frame(night), minn = 50)
#> ibench: 12 officers, 12 tested; flagged at a false discovery rate of 0.05: 2
#>  off   n      rate     bench      ess          z            p            q
#>  o01 252 0.4484127 0.4393818 1738.745  0.2694671 7.875702e-01 9.391220e-01
#>  o02 246 0.4471545 0.4379157 1725.956  0.2727291 7.850615e-01 9.391220e-01
#>  o03 255 0.5294118 0.3417565 2701.527  5.7630501 8.260725e-09 9.912870e-08
#>  o04 266 0.3646617 0.3691144 2733.526 -0.1439985 8.855016e-01 9.391220e-01
#>  o05 247 0.3805668 0.3619805 2740.723  0.5767080 5.641367e-01 9.391220e-01
#>  o06 273 0.3040293 0.3701072 2717.313 -2.2521013 2.431587e-02 9.726348e-02
#>  o07 246 0.3252033 0.3518533 2638.673 -0.8519623 3.942350e-01 9.391220e-01
#>  o08 235 0.3531915 0.3488650 2641.215  0.1330054 8.941891e-01 9.391220e-01
#>  o09 254 0.3425197 0.3449036 2577.094 -0.0763735 9.391220e-01 9.391220e-01
#>  o10 235 0.2936170 0.3549139 2651.348 -1.9692089 4.892911e-02 1.467873e-01
#>  o11 249 0.3694779 0.3490945 2646.200  0.6377712 5.236226e-01 9.391220e-01
#>  o12 242 0.2768595 0.3613048 2687.908 -2.7944882 5.198194e-03 3.118917e-02
#>   flag
#>  FALSE
#>  FALSE
#>   TRUE
#>  FALSE
#>  FALSE
#>  FALSE
#>  FALSE
#>  FALSE
#>  FALSE
#>  FALSE
#>  FALSE
#>   TRUE
ibench(y, off, data.frame(night), minn = 50, dr = TRUE)$tab
#>    off   n      rate     bench      ess           z            p            q
#> 1  o01 252 0.4484127 0.4393818 1738.745  0.27653108 7.821402e-01 9.384964e-01
#> 2  o02 246 0.4471545 0.4379157 1725.956  0.28607426 7.748212e-01 9.384964e-01
#> 3  o03 255 0.5294118 0.3417565 2701.527  5.92119984 3.196012e-09 3.835214e-08
#> 4  o04 266 0.3646617 0.3691144 2733.526 -0.14637442 8.836258e-01 9.384964e-01
#> 5  o05 247 0.3805668 0.3619805 2740.723  0.58858288 5.561411e-01 9.384964e-01
#> 6  o06 273 0.3040293 0.3701072 2717.313 -2.27360614 2.298968e-02 9.195871e-02
#> 7  o07 246 0.3252033 0.3518533 2638.673 -0.86671658 3.860973e-01 9.266336e-01
#> 8  o08 235 0.3531915 0.3488650 2641.215  0.13672962 8.912445e-01 9.384964e-01
#> 9  o09 254 0.3425197 0.3449036 2577.094 -0.07715986 9.384964e-01 9.384964e-01
#> 10 o10 235 0.2936170 0.3549139 2651.348 -2.01670530 4.372628e-02 1.311788e-01
#> 11 o11 249 0.3694779 0.3490945 2646.200  0.65262224 5.139998e-01 9.384964e-01
#> 12 o12 242 0.2768595 0.3613048 2687.908 -2.86006053 4.235602e-03 2.541361e-02
#>     flag
#> 1  FALSE
#> 2  FALSE
#> 3   TRUE
#> 4  FALSE
#> 5  FALSE
#> 6  FALSE
#> 7  FALSE
#> 8  FALSE
#> 9  FALSE
#> 10 FALSE
#> 11 FALSE
#> 12  TRUE
```
