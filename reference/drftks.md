# Two-sample Kolmogorov-Smirnov test (C backend)

Tests whether two numeric samples come from the same distribution, for
example a column of last year's release of a table and the same column
of this year's. It needs no assumption about the shape of either
distribution and reacts to any kind of change: a shift, a change in
spread or a change in shape.

## Usage

``` r
drftks(x, y)
```

## Arguments

- x, y:

  Numeric vectors, the reference and the new sample; logical vectors are
  read as 0/1.

## Value

A named numeric vector of class `drftks`: `statistic` (the KS \\D\\),
`p_value` and `neff` (\\n_xn_y/(n_x+n_y)\\), with the sample sizes in
the attribute `lakner`.

## Details

The statistic is the largest vertical gap between the two empirical
distribution functions, with the asymptotic two-sided p-value \\2\sum_k
(-1)^{k-1} e^{-2k^2t^2}\\, \\t = \sqrt{n\_{\mathrm{eff}}}D\\, computed
for \\t \< 1\\ by the theta-function form, which converges there. It is
accurate for moderate samples and conservative for small ones; for an
exact small-sample p-value use
[`stats::ks.test()`](https://rdrr.io/r/stats/ks.test.html). Ties are
handled by comparing the two functions at each distinct value, so tied
data gives no warning.

Missing values are dropped. Text, a factor or a list, an infinite value,
or no value left in either sample is refused.

## See also

[`drftpsi()`](https://rootcoder007.github.io/crimlab/reference/drftpsi.md),
[`drftchisq()`](https://rootcoder007.github.io/crimlab/reference/drftchisq.md),
[`drfthomog()`](https://rootcoder007.github.io/crimlab/reference/drfthomog.md).

## Examples

``` r
set.seed(1)
ref <- stats::rnorm(200)
drftks(ref, stats::rnorm(200))
#> drftks: statistic =   0.0600, p_value =   0.8643, neff = 100.0000 
#> inputs: nx = 200, ny = 200 
drftks(ref, stats::rnorm(200, mean = 0.8))
#> drftks: statistic = 3.250e-01, p_value = 1.338e-09, neff = 1.000e+02 
#> inputs: nx = 200, ny = 200 
drftks(ref, stats::rnorm(200, sd = 2.5))
#> drftks: statistic = 3.150e-01, p_value = 4.813e-09, neff = 1.000e+02 
#> inputs: nx = 200, ny = 200 
drftks(ref, ref)[["statistic"]]
#> [1] 0
# simulated monthly counts of an offence in two releases
drftks(c(41, 38, 52, 47, 44, 50), c(41, 38, 52, 47, 61, 66))
#> drftks: statistic = 0.3333, p_value = 0.8928, neff = 3.0000 
#> inputs: nx = 6, ny = 6 
```
