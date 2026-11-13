# Fairness metrics of a binary risk tool by group

The confusion matrix of a risk score cut at `cut` within each group, and
the usual rates computed from it: accuracy, sensitivity (true positive
rate), specificity, false positive and false negative rates, positive
and negative predictive values, the share and number predicted positive
and the Matthews correlation, each with its ratio to the base group and,
for the first eight, a Wilson score interval (Wilson 1927). A case is
predicted positive when its score exceeds `cut`. No single metric
defines fairness, and when base rates differ between groups the error
rates and the predictive values cannot all be equal at once
(Chouldechova 2017); the table shows which differ, not which should not.

## Usage

``` r
fair(score, y, grp, cut = 0.5, base = 1, conf = 0.95)
```

## Arguments

- score:

  Risk scores: finite numbers, one per case.

- y:

  Outcomes: 0 or 1 (or `FALSE`/`TRUE`), one per case.

- grp:

  Group labels: numbers, text or a factor, one per case, no missing
  values. Groups are ordered as the factor's levels, or sorted.

- cut:

  The cut: one number.

- base:

  The group the ratios are taken against: its label, or its position
  from 1.

- conf:

  The confidence level of the intervals, strictly between 0 and 1.

## Value

A list of data frames, one row per group with its label in `group`:
`cnt` (`tp`, `fp`, `tn`, `fn`, `n`), `met` (`acc`, `tpr`, `tnr`, `fpr`,
`fnr`, `ppv`, `npv`, `ppos`, `npos`, `mcc`; `NaN` where a rate has an
empty denominator), `rat` (each metric over the base group's), `lo` and
`hi` (the interval ends of the first eight metrics; the last two columns
`NaN`).

## References

Wilson, E. B. (1927). Probable inference, the law of succession, and
statistical inference. *Journal of the American Statistical
Association*, 22(158), 209-212.
[doi:10.1080/01621459.1927.10502953](https://doi.org/10.1080/01621459.1927.10502953)

Chouldechova, A. (2017). Fair prediction with disparate impact: a study
of bias in recidivism prediction instruments. *Big Data*, 5(2), 153-163.
[doi:10.1089/big.2016.0047](https://doi.org/10.1089/big.2016.0047)

## Examples

``` r
set.seed(1)
g <- rep(c("a", "b"), each = 100)
y <- rbinom(200, 1, ifelse(g == "a", 0.3, 0.5))
s <- plogis(rnorm(200, ifelse(y == 1, 0.8, -0.8)))
fair(s, y, g)$met
#>   group acc       tpr       tnr       fpr       fnr       ppv       npv ppos
#> 1     a 0.8 0.8125000 0.7941176 0.2058824 0.1875000 0.6500000 0.9000000 0.40
#> 2     b 0.8 0.8148148 0.7826087 0.2173913 0.1851852 0.8148148 0.7826087 0.54
#>   npos       mcc
#> 1   40 0.5776155
#> 2   54 0.5974235
```
