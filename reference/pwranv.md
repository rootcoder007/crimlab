# Power of the one-way analysis of variance

The power of the F test of equal means in `k` groups of `n` with Cohen's
effect size `f` (the standard deviation of the group means about their
mean, over the common standard deviation), from the noncentral F
distribution; or, given the power wanted, the smallest whole `n` per
group that reaches it.

## Usage

``` r
pwranv(k, n = NULL, f, alpha = 0.05, power = NULL)
```

## Arguments

- k:

  The number of groups: a whole number from 2 to 1000000.

- n:

  `NULL`, or the sizes of each group to try: numbers above 1, at most
  10000000.

- f:

  The effect sizes to try: numbers from 0.

- alpha:

  The level of the test: one number strictly between 0 and 1.

- power:

  `NULL`, or the power wanted: one number above `alpha` and below 1.
  Give one of `n` and `power`.

## Value

A data frame, one row per `n` and `f` (`n` varying fastest), or per `f`
with `power`: `k`, `n`, `f`, `df1`, `df2`, `ncp`, `fcrit` and `power`.

## Details

The statistic has `df1 = k - 1` and `df2 = k (n - 1)` degrees of freedom
and noncentrality `ncp = k n f^2` (Cohen 1988, ch. 8); the power is the
chance that it passes `fcrit`, the upper `alpha` quantile of the central
F. The noncentral upper tail is the Poisson mixture of beta upper tails
(Johnson, Kotz and Balakrishnan 1995, ch. 30). With `power`, `n` is
searched over the whole numbers from 2 to 10000000; where no such `n`
reaches it (`f = 0`, or too small an effect) the row is `NA` from `df2`
on. Cohen's conventions are 0.1, 0.25 and 0.4 for a small, medium and
large effect.

## References

Cohen, J. (1988). *Statistical Power Analysis for the Behavioral
Sciences*, 2nd ed. Lawrence Erlbaum.
[doi:10.4324/9780203771587](https://doi.org/10.4324/9780203771587)

Johnson, N. L., Kotz, S. and Balakrishnan, N. (1995). *Continuous
Univariate Distributions*, vol. 2, 2nd ed. Wiley.

## Examples

``` r
pwranv(3, n = c(20, 40, 60), f = 0.25)
#>   k  n    f df1 df2   ncp    fcrit     power
#> 1 3 20 0.25   2  57  3.75 3.158843 0.3744311
#> 2 3 40 0.25   2 117  7.50 3.073763 0.6756745
#> 3 3 60 0.25   2 177 11.25 3.047012 0.8546381
pwranv(4, f = c(0.1, 0.25, 0.4), power = 0.8)
#>   k   n    f df1  df2   ncp    fcrit     power
#> 1 4 274 0.10   3 1092 10.96 2.613053 0.8007324
#> 2 4  45 0.25   3  176 11.25 2.655939 0.8039869
#> 3 4  19 0.40   3   72 12.16 2.731807 0.8234006
```
