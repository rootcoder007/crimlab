# Lengths of stay: their distribution and its fit

The lengths of stay of a set of placements (or of people): the mean with
its t interval, the median and quartiles, the longest; the share of
stays in each band and at or below each threshold (15 days, the Mandela
Rules' limit, say), with the share of each band carrying a flag (no
review on record); and the fit of lognormal, exponential and Pareto laws
by maximum likelihood, compared by AIC: a heavy tail (Pareto or
lognormal beating the exponential) separates a few long stays from many
short ones. A stay still open is a lower bound on its length; see
[`losbnd()`](https://rootcoder007.github.io/crimlab/reference/losbnd.md)
for the mean it bounds.

## Usage

``` r
lostay(days, w = NULL, breaks = NULL, thr = NULL, flag = NULL, level = 0.95)
```

## Arguments

- days:

  The stays: numbers from 0.

- w:

  `NULL`, or how many stays each value stands for: whole numbers from 0
  with a positive total.

- breaks:

  `NULL`, or band edges, increasing; a band is `(lo, hi]`.

- thr:

  `NULL`, or thresholds: numbers from 0.

- flag:

  `NULL`, or `TRUE`/`FALSE` for each stay.

- level:

  The level of the interval of the mean: one number between 0 and 1.

## Value

An object of class `lostay`, a list: `stay`, a named vector (`n`,
`days`, `mean`, `sd`, `median`, `q1`, `q3`, `max`, `se`, `lo`, `hi`);
`bands` (with `breaks`: `lo`, `hi`, `n`, `share` of all stays, with
`flag` `flagged` the share of the band's stays); `thr` (`thr`, `atmost`,
`above`); `fit` (`model`, `aic`, `best`; the fits use the stays above
0); and `par` (`meanlog`, `sdlog`, `alpha`).

## References

Clauset, A., Shalizi, C. R. and Newman, M. E. J. (2009). Power-law
distributions in empirical data. *SIAM Review*, 51(4), 661-703.
[doi:10.1137/070710111](https://doi.org/10.1137/070710111)

## Examples

``` r
set.seed(3)
s <- round(rlnorm(400, 2.3, 1.1)) + 1
r <- lostay(s, breaks = c(0, 5, 15, 31, 61, Inf), thr = 15)
r$bands
#>   lo  hi   n  share
#> 1  0   5  93 0.2325
#> 2  5  15 157 0.3925
#> 3 15  31  84 0.2100
#> 4 31  61  43 0.1075
#> 5 61 Inf  23 0.0575
r$fit
#>         model      aic  best
#> 1   lognormal 3105.502  TRUE
#> 2 exponential 3228.711 FALSE
#> 3      pareto 3249.590 FALSE
```
