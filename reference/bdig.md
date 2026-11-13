# Benford digit screens

`bdig()` takes the first `nd` significant digits of each value. With
`dp = -1` they are the digits of the shortest decimal that gives back
the double, decided exactly, so 0.29 gives 29 where truncating `x` times
a power of ten gives 28; with `dp` from 0 to 22 the data are taken as
recorded to `dp` decimals. `bnft()` tests the digits against Benford's
law, `P(d) = log10(1 + 1 / d)`, by the chi-square statistic with its
exact upper tail, a Monte Carlo p-value from samples of Benford digits
of the same size, and Nigrini's mean absolute deviation. A departure
from Benford's law is a reason to look, not a finding of fabrication:
many honest data sets (bounded amounts, assigned numbers) do not follow
it.

## Usage

``` r
bdig(x, nd = 1, dp = -1, sgn = "pos")

bnft(dig, nd = 1, nsim = 0, seed = NULL)
```

## Arguments

- x:

  Numbers. Zeros, non-finite values (including `NA`), values of the
  excluded sign and values outside the exact range are set aside and
  counted.

- nd:

  The number of leading digits: from 1 to 9 for `bdig()`, 1 to 6 for
  `bnft()`.

- dp:

  `-1`, or the number of decimals the data were recorded to.

- sgn:

  Which values to keep: `"pos"`, `"neg"` or `"both"` (by absolute
  value).

- dig:

  Digits as `bdig()` returns them (0 is skipped).

- nsim:

  The number of Monte Carlo samples: a whole number from 0.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

`bdig()`: a list of `dig` (an integer per value, 0 when set aside) and
`cnt` (zeros, non-finite, excluded sign, out of range). `bnft()`: a list
of `obs` (a data frame of `digit`, observed `count` and Benford
`exp`ected count) and `stat` (`m` digits used, `chisq`, `df`, `p`, `pmc`
(Monte Carlo, `NaN` without samples) and `mad`).

## References

Benford, F. (1938). The law of anomalous numbers. *Proceedings of the
American Philosophical Society*, 78(4), 551-572.

Nigrini, M. J. (2012). *Benford's Law: Applications for Forensic
Accounting, Auditing, and Fraud Detection*. Wiley.
[doi:10.1002/9781119203094](https://doi.org/10.1002/9781119203094)

## Examples

``` r
set.seed(1)
x <- exp(runif(500, 0, 10))
d <- bdig(x)$dig
bnft(d, nsim = 199, seed = 1)$stat
#>            m        chisq           df            p          pmc          mad 
#> 500.00000000  15.41054423   8.00000000   0.05163739   0.05000000   0.01580472 
```
