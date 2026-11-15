# Chi-square test of a categorical column against a known distribution

Tests whether the category counts of a new sample agree with a reference
distribution taken as known (proportions fixed by a specification, or a
published breakdown): Pearson's chi-square goodness-of-fit test. When
the reference is itself a sample, use
[`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md),
which carries the reference's own sampling error.

## Usage

``` r
drftchisq(obs, expd)
```

## Arguments

- obs:

  Named counts of the new sample, or a factor or text to be tabulated.

- expd:

  Named reference counts or proportions, or a factor or text to be
  tabulated; rescaled to the total of `obs`.

## Value

A named numeric vector of class `drftchisq`: `statistic`, `df`,
`p_value` and `logpvalue` (the log of the upper tail, which keeps its
information where `p_value` underflows), with the inputs in the
attribute `lakner`.

## Details

Categories present in one argument and not the other are aligned by
name, so a vanished category registers as a shortfall. With \\O_j\\ the
observed counts and \\E_j = N e_j / \sum e\\ the reference scaled to the
observed total, \$\$X^2 = \sum_j \frac{(O_j - E_j)^2}{E_j}\$\$ over the
categories with \\E_j \> 0\\, on \\k - 1\\ degrees of freedom; the
p-value is the upper tail, computed directly. A category that occurs but
which the reference gives probability zero contradicts it outright:
`statistic` is then `Inf` and `p_value` 0. With a single category there
is nothing to test: `df` is 0 and the p-values `NA`. The approximation
is poor below expected counts of about 5.

## See also

[`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md),
[`stats::chisq.test()`](https://rdrr.io/r/stats/chisq.test.html).

## Examples

``` r
ref <- c(a = 50, b = 30, c = 20)
drftchisq(c(a = 100, b = 60, c = 40), ref)
#> drftchisq: statistic = 0, df = 2, p_value = 1, logpvalue = 0 
#> inputs: n = 200, k = 3 
drftchisq(c(a = 40, b = 60, c = 100), ref)
#> drftchisq: statistic =  1.26e+02, df =  2.00e+00, p_value =  4.36e-28, logpvalue = -6.30e+01 
#> inputs: n = 200, k = 3 
drftchisq(c("a", "a", "b", "b"), c("a", "b"))
#> drftchisq: statistic = 0, df = 1, p_value = 1, logpvalue = 0 
#> inputs: n = 4, k = 2 
drftchisq(c("a", "a", "b", "d"), c("a", "a", "b", "b"))
#> drftchisq: statistic =  Inf, df =    1, p_value =    0, logpvalue = -Inf 
#> inputs: n = 4, k = 3 
```
