# Bayes outbreak detection

The Bayes detector of `surveillance::algo.bayes()` (Riebler 2004; Hoehle
2007) for a series of counts: at each time `t` the reference counts are
the `w` before it (`acty = TRUE`) and, in each of the `b` years before,
the `2 w + 1` counts around the same time of year (`freq` counts per
year); missing counts are left out. With `s` their sum and `n` their
number, the posterior predictive distribution of the count at `t` under
a Poisson model with a flat prior is negative binomial of size `s + 0.5`
and probability `n / (n + 1)`; its `1 - alpha` quantile is the upper
bound, and a count above it raises an alarm. The same series may be
calls for service, emergency visits or reported crimes of one kind.

## Usage

``` r
bayout(x, freq = 52, b = 0, w = 6, acty = TRUE, alpha = 0.05, at = NULL)
```

## Arguments

- x:

  The counts: numbers from 0, `NA` for a missing count.

- freq:

  The counts per year: one whole number from 1 (52 for weeks).

- b:

  The number of past years: one whole number from 0.

- w:

  The half-width of the window: one whole number from 0.

- acty:

  `TRUE` to use the `w` counts just before `t`.

- alpha:

  The level: one number strictly between 0 and 1.

- at:

  `NULL` for every time from `b freq + w + 1` on, or the times to
  evaluate (from 1).

## Value

A data frame with one row per time: `t`, `ub` (the upper bound) and
`alarm` (`NA` where the count is missing).

## References

Hoehle, M. (2007). surveillance: an R package for the monitoring of
infectious diseases. *Computational Statistics*, 22(4), 571-582.
[doi:10.1007/s00180-007-0074-8](https://doi.org/10.1007/s00180-007-0074-8)

Riebler, A. (2004). *Empirischer Vergleich von statistischen Methoden
zur Ausbruchserkennung bei Surveillance Daten*. Bachelor's thesis,
University of Munich.

## See also

Other spatial epidemiology:
[`bufexp()`](https://rootcoder007.github.io/crimlab/reference/bufexp.md),
[`bym2sc()`](https://rootcoder007.github.io/crimlab/reference/bym2sc.md),
[`ecolrg()`](https://rootcoder007.github.io/crimlab/reference/ecolrg.md),
[`evfilt()`](https://rootcoder007.github.io/crimlab/reference/evfilt.md),
[`kerexp()`](https://rootcoder007.github.io/crimlab/reference/kerexp.md),
[`lerprc()`](https://rootcoder007.github.io/crimlab/reference/lerprc.md),
[`sfrail()`](https://rootcoder007.github.io/crimlab/reference/sfrail.md),
[`spbias()`](https://rootcoder007.github.io/crimlab/reference/spbias.md),
[`spsir()`](https://rootcoder007.github.io/crimlab/reference/spsir.md)

## Examples

``` r
bayout(c(3, 5, 2, 4, 6, 3, 4, 12), w = 6, at = 8)
#>   t ub alarm
#> 1 8  8  TRUE
```
