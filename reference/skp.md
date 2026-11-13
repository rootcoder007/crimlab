# Scott-Knott grouping of measurements

Splits sorted measurements into groups that differ (Scott and Knott
1974) when the variance of one measurement is known, as in the forensic
glass version of Triggs, Curran, Buckleton and Walsh (1997): a segment
is split at the cut with the largest between-groups sum of squares when
that is significant at `alpha`, judged against `nsim` simulated samples
of independent normal values with the known variance (a test of exact
size under that null).

## Usage

``` r
skp(x, s2, alpha = 0.05, nsim = 999, seed = NULL)
```

## Arguments

- x:

  Measurements: finite numbers (refractive indices, say).

- s2:

  The variance of one measurement, on the scale of `x`: positive.

- alpha:

  The level of each split's test, strictly between 0 and 1.

- nsim:

  The number of simulated samples per test: a whole number from 1.

- seed:

  `NULL` to draw the simulations' seed from R's generator (so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes them), or one
  whole number from 0.

## Value

An integer vector of group labels, one per measurement in input order,
numbered from 1 in increasing order of the values.

## References

Scott, A. J. and Knott, M. (1974). A cluster analysis method for
grouping means in the analysis of variance. *Biometrics*, 30(3),
507-512. [doi:10.2307/2529204](https://doi.org/10.2307/2529204)

Triggs, C. M., Curran, J. M., Buckleton, J. S. and Walsh, K. A. J.
(1997). The grouping problem in forensic glass analysis: a divisive
approach. *Forensic Science International*, 85(1), 1-14.
[doi:10.1016/S0379-0738(96)02037-3](https://doi.org/10.1016/S0379-0738%2896%2902037-3)

## Examples

``` r
set.seed(1)
x <- c(rnorm(5, 1.5180, 4e-5), rnorm(5, 1.5195, 4e-5))
skp(x, s2 = 4e-5^2, nsim = 199, seed = 1)
#>  [1] 1 1 1 1 1 2 2 2 2 2
```
