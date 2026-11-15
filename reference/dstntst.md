# Do distinct groups differ? A heterogeneity test before pooling

Tests whether the estimates of one indicator for distinct groups (First
Nations, Metis and Inuit; Aboriginal and Torres Strait Islander people;
any groups with published intervals) differ more than their intervals
allow, and says how much of their spread is real: the evidence for
reporting each group distinctly rather than as one pooled figure.

## Usage

``` r
dstntst(
  data,
  grp,
  est,
  lower,
  upper,
  qual = NULL,
  conflvl = 0.95,
  pooled = NULL
)

# S3 method for class 'dstntst'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per group.

- grp, est, lower, upper:

  Column names.

- qual:

  Optional column of quality flags; suppressed rows are left out.

- conflvl:

  Confidence level of the published intervals, and the level of the
  verdict printed (default 0.95).

- pooled:

  Optional labels of rows that pool other groups.

- x:

  A result of `dstntst()`.

- ...:

  Not used.

## Value

A list of class `dstntst`: `groups` (the estimates used, with their
standard errors), `pooled` (the inverse-variance mean) and `poose`, `q`,
`df`, `p_value`, `i2`, `leftout` (the rows not used) and `conflvl`.
[`print()`](https://rdrr.io/r/base/print.html) shows it.

## Details

With estimates \\y_g\\, standard errors \\s_g\\ from the published
intervals and weights \\w_g = 1/s_g^2\\, Cochran's \\Q = \sum w_g (y_g -
\bar y_w)^2\\ has a chi-squared distribution on \\G - 1\\ degrees of
freedom when the groups share one value, and \\I^2 = \max(0, (Q -
(G-1))/Q)\\ is the share of the spread that is real rather than sampling
error (Higgins and Thompson 2002). A pooled figure is a weighted mean of
different populations whenever \\I^2\\ is large; a population-weighted
total also depends on how many people are in each group, so it is mostly
the largest group's figure. Rows that pool groups are left out
(`pooled`, or labels starting with "Total", "Indigenous identity",
"Aboriginal identity", "Total Indigenous" or "Indigenous population");
every other row is included, so drop a reference group yourself, and
never test overlapping groups (`"alone"` beside “alone or in
combination”) together: they are not independent. Rows with a quality
flag `"F"`, `"x"`, `".."` or `"..."` (suppressed) are left out; other
flags than those and `"E"` are treated as usable, with a warning. An
estimate outside its own interval is not compared, with a warning. With
two or three groups the test has little power and \\I^2\\ is imprecise;
a large \\p\\ is not evidence that the groups are the same.

## References

Higgins, J. P. T. and Thompson, S. G. (2002). Quantifying heterogeneity
in a meta-analysis. Statistics in Medicine 21(11), 1539-1558.
[doi:10.1002/sim.1186](https://doi.org/10.1002/sim.1186)

## Examples

``` r
# your own table: made-up numbers, with the pooled row named
own <- data.frame(
  group = c("Aboriginal", "Torres Strait Islander", "Both"),
  est = c(24.0, 31.0, 24.6),
  lo = c(22.1, 26.2, 22.8),
  hi = c(25.9, 35.8, 26.4)
)
dstntst(own, "group", "est", "lo", "hi", pooled = "Both")
#> Distinctions test: 2 groups
#>   Cochran's Q 7.06 on 1 df, p = 0.00787
#>   I-squared: 86% of the variation is beyond chance
#>   Pooled: 24.95 (se 0.901, inverse-variance mean)
#>   Left out: Both
#>                   group estimate        se
#>              Aboriginal       24 0.9694056
#>  Torres Strait Islander       31 2.4490246

# simulated victimisation rates of three police divisions
div <- data.frame(
  division = c("North", "Central", "South"),
  rate = c(4.1, 5.3, 4.4), lo = c(3.6, 4.7, 3.9), hi = c(4.6, 5.9, 4.9)
)
dstntst(div, "division", "rate", "lo", "hi")
#> Distinctions test: 3 groups
#>   Cochran's Q 9.42 on 2 df, p = 0.00899
#>   I-squared: 79% of the variation is beyond chance
#>   Pooled: 4.52 (se 0.155, inverse-variance mean)
#>    group estimate        se
#>    North      4.1 0.2551067
#>  Central      5.3 0.3061281
#>    South      4.4 0.2551067
```
