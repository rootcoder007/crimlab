# Differences between groups in a risk tool's decisions

How the share of cases a tool flags (predicts positive) differs between
groups, as a difference and as a ratio against a base group, each with
an interval; with the outcomes, how its true and false positive rates
differ; and the bias amplification score, the largest parity gap times
the Gini coefficient of the groups' shares.

## Usage

``` r
fairdf(
  pred,
  grp = NULL,
  y = NULL,
  base = NULL,
  conf = 0.95,
  nboot = 0,
  seed = NULL
)
```

## Arguments

- pred:

  The tool's decisions: 1 (or `TRUE`) for a case it flags, 0 (or
  `FALSE`) otherwise, one per case. A `clpre` list from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too: its elements fill the arguments by name.

- grp:

  Group labels: numbers, text or a factor, one per case, no missing
  values. Groups are ordered as the factor's levels, or sorted.

- y:

  `NULL`, or the outcomes: 0 or 1 (or `FALSE`/`TRUE`), one per case.

- base:

  The base group: `NULL` for the group with the largest share flagged
  (the first in order on ties), its label, or its position from 1.

- conf:

  The confidence level of the intervals, strictly between 0 and 1.

- nboot:

  The bootstrap replicates for the intervals of the Gini coefficient and
  the score: a whole number from 0 (none).

- seed:

  `NULL` to draw the seed of the bootstrap (or the permutations) from
  R's generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it, or one
  whole number from 0.

## Value

A list: `tab`, a data frame with one row per group (in the order of
`grp`): `group`, `n`, `npos` (cases flagged), `ppos` (their share) with
`lo` and `hi`, `dp` with `dplo` and `dphi`, `di` with `dilo` and `dihi`,
and with `y` also `tpr`, `fpr`, `dtpr` (`dtprlo`, `dtprhi`), `dfpr`
(`dfprlo`, `dfprhi`), `aod` (`aodlo`, `aodhi`) and `eo`; the base
group's own intervals are `NaN`. `sum`, a named vector: `di` (the
smallest finite ratio of another group), `dp` (the gap of largest size,
signed), `gini`, `bas`, `ginilo`, `ginihi`, `baslo`, `bashi` (`NaN` when
`nboot` is 0) and with `y` also `aod` and `eo` (those of largest size,
signed). `base`, the base group's label.

## Details

[`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md)
reports each group's rates and their ratios; `fairdf()` reports the
comparisons the fairness literature names. The parity gap `dp` is a
group's share flagged minus the base group's (demographic parity holds
when every gap is 0), with Newcombe's (1998) hybrid score interval. The
ratio `di` is the disparate impact ratio, with the log interval of Katz
et al. (1978); under the four-fifths rule a ratio below 0.8 against the
group with the highest selection rate is the usual sign of adverse
impact, so by default that group is the base. With the outcomes `y`,
`dtpr` and `dfpr` are the gaps in true and false positive rates (both 0
is equalized odds), `aod` their mean (the average odds difference) with
the interval of Zou and Donner (2008) built from the two Newcombe
intervals, and `eo` the larger of the two gaps in size. The bias
amplification score of Barman and Barman (2026) is large only when a
sizeable gap comes with unequal shares across all the groups; its
bootstrap interval, and the Gini coefficient's, are percentile intervals
from redrawing each group's decisions with replacement (Efron and
Tibshirani 1993), the base group held.

These numbers describe what the tool did to each group. They do not say
why: a gap in shares flagged can come from a gap in the outcomes the
tool was built to predict, from how those outcomes were recorded, or
from the tool; and when the groups' base rates differ, equal error rates
and equal predictive values cannot all hold at once (Chouldechova 2017,
see
[`impfpr()`](https://rootcoder007.github.io/crimlab/reference/impfpr.md)).
A ratio's interval is `NaN` when either share is 0, and a ratio against
a base group with no case flagged is `Inf` (or `NaN` for 0 over 0).

## References

Newcombe, R. G. (1998). Interval estimation for the difference between
independent proportions: comparison of eleven methods. *Statistics in
Medicine*, 17(8), 873-890.
[doi:10.1002/(SICI)1097-0258(19980430)17:8\<873::AID-SIM779\>3.0.CO;2-I](https://doi.org/10.1002/%28SICI%291097-0258%2819980430%2917%3A8%3C873%3A%3AAID-SIM779%3E3.0.CO%3B2-I)

Katz, D., Baptista, J., Azen, S. P. and Pike, M. C. (1978). Obtaining
confidence intervals for the risk ratio in cohort studies. *Biometrics*,
34(3), 469-474. [doi:10.2307/2530610](https://doi.org/10.2307/2530610)

Zou, G. Y. and Donner, A. (2008). Construction of confidence limits
about effect measures: a general approach. *Statistics in Medicine*,
27(10), 1693-1702.
[doi:10.1002/sim.3095](https://doi.org/10.1002/sim.3095)

Barman, P. K. and Barman, P. K. (2026). Spatial disparities in simulated
patrol detection: a reproducible generative-model audit across Baltimore
and Chicago.
[doi:10.48550/arXiv.2603.18987](https://doi.org/10.48550/arXiv.2603.18987)

Efron, B. and Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman and Hall.
[doi:10.1007/978-1-4899-4541-9](https://doi.org/10.1007/978-1-4899-4541-9)

## Examples

``` r
set.seed(1)
g <- rep(c("a", "b"), each = 200)
y <- rbinom(400, 1, 0.3)
p <- rbinom(400, 1, ifelse(g == "a", 0.4, 0.25))
r <- fairdf(p, g, y, nboot = 199, seed = 1)
r$tab[, c("group", "ppos", "dp", "dplo", "dphi", "di", "dilo", "dihi")]
#>   group  ppos     dp       dplo        dphi        di      dilo      dihi
#> 1     a 0.435  0.000        NaN         NaN 1.0000000       NaN       NaN
#> 2     b 0.260 -0.175 -0.2640056 -0.08182849 0.5977011 0.4507581 0.7925463
r$sum
#>           di           dp         gini          bas       ginilo       ginihi 
#>  0.597701149 -0.175000000  0.125899281 -0.022032374  0.059494196  0.194227687 
#>        baslo        bashi          aod           eo 
#> -0.054404735 -0.004745647 -0.159044465 -0.202023149 
```
