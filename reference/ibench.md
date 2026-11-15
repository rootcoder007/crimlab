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
  (The paper adds a regression adjustment to the weights, a doubly
  robust estimate, and estimates the local false discovery rate; the
  plain weights and the step-up rate are used here.)

## Usage

``` r
ibench(y, off, x = NULL, minn = 20, fdr = 0.05)
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

## Value

An object of class `ibench`, a list: `tab`, one row per officer (`off`,
`n` stops, `rate`, the officer's share of the outcome, `bench`, the
benchmark, `ess`, its effective number of stops, `z`, `p`, `q` and
`flag`); `fdr`; `minn`; and `nx`, the number of columns of `x`.

## Details

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
```
