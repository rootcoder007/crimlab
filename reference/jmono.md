# Do judges' leniency rankings hold? Monotonicity checks

The decision rate of each judge (with a chi-square test that the rate is
the same for all, and its range), and two checks of the monotonicity a
judge-leniency design rests on.

## Usage

``` r
jmono(
  judge,
  d,
  y = NULL,
  grp = NULL,
  w = NULL,
  lo = NULL,
  hi = NULL,
  minn = 30,
  adj = "holm",
  alpha = 0.05
)
```

## Arguments

- judge:

  The judge of each case: labels.

- d:

  The decision of each case: 0 or 1 (or `FALSE`, `TRUE`).

- y:

  `NULL`, or the outcome of each case: numbers.

- grp:

  `NULL`, or the group of each case: labels.

- w:

  `NULL`, or a weight from 0 for each case (for the rates and outcome
  means; the tests use the counts).

- lo, hi:

  `NULL` (the smallest and largest outcome) or the bounds of the
  outcome.

- minn:

  The fewest cases a judge needs in a group to be compared there: one
  number from 0.

- adj:

  The adjustment of the p-values: `"none"`, `"bonferroni"`, `"holm"`,
  `"hochberg"`, `"BH"` or `"BY"`.

- alpha:

  The level at which an adjusted p-value marks a reversal as
  significant.

## Value

An object of class `jmono`, a list: `judges` (`judge`, `n`, `rate`, with
`y` `outcome`; most lenient first); `homog` (`chisq`, `df`, `p`, `min`,
`max`); `slope`, with `y`, one row per pair (`j` the more lenient, `k`,
`dp`, `dy`, `bound`, `viol`, `late`); `sub`, with `grp`, a list of
`pairs` (`group`, `stricter`, `lenient`, their overall rates `rates`,
`ratel`, their group rates `grates`, `gratel`, their cases `ns`, `nl`,
`z`, `p`, `padj`, `reversed`, `signif`), `summary` (`pairs`, `reversed`,
`share`, `signif`, `refuted` the pairs significantly reversed, `judges`
in them) and `slopes` (`group`, `judges`, `cases`, `slope`); and
`alpha`. The pairs grow as the square of the judges.

## Details

- With outcomes `y` in `[lo, hi]`: if judges are nested (a stricter
  judge detains everyone a more lenient one detains), the outcome means
  of two judges differ by at most `(hi - lo)` times the difference in
  their rates. Every pair is checked; a violation refutes the nesting of
  that pair. The ratio `late` is the effect on the marginal cases
  between them.

- With groups of cases `grp`: nested judges are ordered the same way in
  every group. For each pair of judges with at least `minn` cases in a
  group, the stricter judge overall should not be more lenient there;
  `z` is the pooled two-proportion statistic of the lenient judge's
  group rate over the stricter's, its one-sided p-value is adjusted over
  all comparisons, and `slopes` regresses each group's rates on the
  overall ones (Frandsen, Lefgren and Leslie 2023).

See `proof("jmono")`.

## References

Frandsen, B. R., Lefgren, L. J. and Leslie, E. C. (2023). Judging judge
fixed effects. *American Economic Review*, 113(1), 253-277.
[doi:10.1257/aer.20201860](https://doi.org/10.1257/aer.20201860)

## Examples

``` r
set.seed(9)
j <- sample(paste0("J", 1:6), 1200, TRUE)
len <- c(J1 = 0.2, J2 = 0.3, J3 = 0.35, J4 = 0.4, J5 = 0.5, J6 = 0.6)
d <- rbinom(1200, 1, len[j])
g <- sample(c("adult", "youth"), 1200, TRUE)
r <- jmono(j, d, grp = g)
r$sub$summary
#>       pairs    reversed       share      signif     refuted      judges 
#> 30.00000000  1.00000000  0.03333333  0.00000000  0.00000000  0.00000000 
```
