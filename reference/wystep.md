# Westfall-Young stepdown p-values

Adjusted p-values for several hypotheses tested together, from the joint
null distribution of their statistics (permutations or bootstrap draws),
holding the family-wise error rate while using the correlation between
the tests.

## Usage

``` r
wystep(obs, null)
```

## Arguments

- obs:

  The observed statistics: finite numbers, one per hypothesis.

- null:

  The statistics under the null: a numeric matrix, one row per draw and
  one column per hypothesis.

## Value

A data frame, one row per hypothesis in the order of `obs`: `stat`, `p`
(the share of draws of that statistic at least as large in absolute
value) and `padj` (the stepdown p-value).

## Details

The statistics are ordered by size, `|t_(1)| >= ... >= |t_(K)|`. For
each null draw `b`, `q_bj = max_{r >= j} |t*_b(r)|`; the adjusted
p-value of the `j`-th is the share of draws with `q_bj >= |t_(j)|`, made
monotone: `p~_(j) = max(p~_(j-1), p_(j))` (Westfall and Young 1993,
algorithm 2.8). Draws with a missing or infinite statistic are left out.

## References

Westfall, P. H. and Young, S. S. (1993). *Resampling-Based Multiple
Testing*. Wiley, New York. ISBN 978-0-471-55761-6.

## Examples

``` r
set.seed(1)
null <- matrix(rnorm(3000), ncol = 3)
wystep(c(2.5, 0.4, -2.1), null)
#>   stat     p  padj
#> 1  2.5 0.012 0.044
#> 2  0.4 0.694 0.694
#> 3 -2.1 0.040 0.081
```
