# Chi-square test of homogeneity of two categorical samples

Tests whether two samples of a categorical variable (offence type,
region, race in last year's and this year's release) were drawn from the
same distribution, by Pearson's chi-square on the 2-by-k table of their
counts.

## Usage

``` r
drfthomog(x, y, seed = 1L)
```

## Arguments

- x, y:

  Factors or text (or named count vectors): the reference and the new
  sample.

- seed:

  The seed of the Monte Carlo branch, one whole number from -2147483647
  to 2147483647.

## Value

A named numeric vector of class `drfthomog`: `statistic`, `df`,
`p_value` and `logpvalue`, with the attribute `method` and the inputs in
the attribute `lakner`.

## Details

With pooled shares \\\hat\pi_j = (x_j + y_j)/(n_x + n_y)\\, \$\$X^2 =
\sum_j \frac{(x_j - n_x\hat\pi_j)^2}{n_x\hat\pi_j} + \frac{(y_j -
n_y\hat\pi_j)^2}{n_y\hat\pi_j}\$\$ on \\k - 1\\ degrees of freedom.
Three branches, named in the attribute `method`:
`"chi-square approximation"` when every expected count is at least 5;
`"Monte Carlo permutation (small expected counts)"` otherwise, the
p-value from
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html) with
`simulate.p.value = TRUE` and `B = 10000` tables with the same margins,
drawn under `seed` with the caller's random stream restored (`p_value`
is then at least 1/(B+1)); and `"inapplicable: ..."` for a single
category or fewer than two observations per category, where no test has
power and `p_value` is `NA` rather than 1.

[`drftchisq()`](https://rootcoder007.github.io/crimlab/reference/drftchisq.md)
takes the reference as known; when it is itself a finite sample, that
ignores its sampling error and reports drift too readily. This test
estimates the shared distribution from the pooled margins and carries
the uncertainty of both samples, the right test for two releases of a
table. Categories present in only one sample are given zero counts in
the other.

## See also

[`drftchisq()`](https://rootcoder007.github.io/crimlab/reference/drftchisq.md),
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html).

## Examples

``` r
set.seed(1)
a <- sample(c("x", "y", "z"), 300, TRUE)
b <- sample(c("x", "y", "z"), 300, TRUE)
drfthomog(a, b)
#> drfthomog: statistic =  0.5150, df =  2.0000, p_value =  0.7730, logpvalue = -0.2575 
#> inputs: nx = 300, ny = 300, k = 3 
#> method: chi-square approximation 
drfthomog(a, sample(c("x", "y", "z"), 300, TRUE, prob = c(0.7, 0.2, 0.1)))
#> drfthomog: statistic =  6.978e+01, df =  2.000e+00, p_value =  7.030e-16, logpvalue = -3.489e+01 
#> inputs: nx = 300, ny = 300, k = 3 
#> method: chi-square approximation 
drfthomog(c(a = 30, b = 20), c(a = 25, b = 25))
#> drfthomog: statistic =  1.0101, df =  1.0000, p_value =  0.3149, logpvalue = -1.1556 
#> inputs: nx = 50, ny = 50, k = 2 
#> method: chi-square approximation 
attr(drfthomog(c(a = 3, b = 1, c = 2), c(a = 1, b = 4, c = 1)), "method")
#> [1] "Monte Carlo permutation (small expected counts)"
try(drfthomog(c(1, 2), c(a = 3)))
#> Error : drfthomog: a numeric `x` must be named.
```
