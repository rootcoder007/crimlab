# Compare published group estimates with a reference group

For published estimates with their confidence limits (percentages,
rates, years of life expectancy: any scale, with the limits on the same
scale), the ratio and the difference of each group to a reference group,
with confidence intervals worked from the published ones: First Nations,
Metis and Inuit against non-Indigenous people, Aboriginal and Torres
Strait Islander people against non-Indigenous Australians, Maori against
non-Maori, or any groups with published intervals.

## Usage

``` r
indgap(
  data,
  grp,
  est,
  lower,
  upper,
  ref = "Non-Indigenous",
  qual = NULL,
  conflvl = 0.95,
  pooled = NULL
)
```

## Arguments

- data:

  A data frame, one row per group (filtered to one period, age group,
  sex and place: a group appearing twice is refused).

- grp, est, lower, upper:

  Column names: the group label, the estimate, and its published
  confidence limits, on one scale.

- ref:

  The label of the reference group (default `"Non-Indigenous"`).

- qual:

  Optional column of quality flags in the package's vocabulary (`"E"`,
  `"F"`, `"x"`, `".."`, `"..."`; see
  [`qltyrcd`](https://rootcoder007.github.io/crimlab/reference/qltyrcd.md)
  for other agencies' flags).

- conflvl:

  The confidence level of the published intervals and of the ones
  returned (default 0.95).

- pooled:

  Optional labels of the rows that pool other groups (for example
  `"Aboriginal and Torres Strait Islander"` beside its parts, or a
  “Total Indigenous” row); by default the labels described in Details.
  Each label must be a group of `data`: one that matches no row (a typo)
  is refused rather than ignored.

## Value

A data frame of class `indgap`, one row per group but the reference. The
columns are `group`, `estimate`, `refest`, `ratio`, `ratiolo`,
`ratioup`, `difference`, `diflo`, `difup`, `usewithcau` and `poogrp`;
attributes `reference` and `conflvl`.

## Details

Each published interval gives a standard error, \\s = (u - l) / (2 z)\\,
on the scale of the estimate for the difference and on the log scale,
\\\sigma = (\log u - \log l) / (2 z)\\, for the ratio, with \\z =
\Phi^{-1}(1 - \alpha/2)\\. This takes the interval as symmetric on that
scale; for asymmetric (gamma, logit, Wilson) limits it averages the two
sides. With a standard error, a coefficient of variation or a margin of
error instead of limits, make limits with
[`limfrom`](https://rootcoder007.github.io/crimlab/reference/limfrom.md)
first. The groups are treated as independent samples, so the standard
error of a difference is \\\sqrt{s_g^2 + s_r^2}\\ and of a log ratio
\\\sqrt{\sigma_g^2 + \sigma_r^2}\\. Disjoint domains of one survey are
nearly independent but not exactly (they share the design's strata and
clusters): the intervals are approximations, and the exact ones need the
survey's replicate weights (for Statistics Canada's surveys, in its
Research Data Centres). The reference must not contain the groups:
compare with non-Indigenous people, not with a total, whose covariance
with each group is large and positive (the function warns when the
reference looks like a total). A lower limit of 0 (a limit cut at zero)
leaves the ratio's interval `NA`.

Distinctions-based: each group is compared on its own. A row that pools
groups is kept but marked `poogrp = TRUE`, since a pooled figure can
hide what differs between the groups it pools: name such rows in
`pooled` (by default labels beginning `"Total"`, “Indigenous identity”,
“Aboriginal identity” or “Indigenous population” are taken as pooled). A
row flagged `E` (in the group or the reference) is compared and marked;
one flagged `F`, `x`, `..` or `...`, with no estimate, or whose
published interval does not contain its own estimate, is not compared.
Other agencies' flags:
[`qltyrcd`](https://rootcoder007.github.io/crimlab/reference/qltyrcd.md).

## References

Statistics Canada (2025). Health indicators for First Nations adults
living off reserve, Metis and Inuit adults, 2015 to 2022. The Daily,
2025-02-13.  
www150.statcan.gc.ca/n1/daily-quotidien/250213/dq250213f-eng.htm

Carroll SR, Garba I, Figueroa-Rodriguez OL, et al. (2020). The CARE
Principles for Indigenous Data Governance. Data Science Journal 19:43.
[doi:10.5334/dsj-2020-043](https://doi.org/10.5334/dsj-2020-043)

## See also

[`limfrom`](https://rootcoder007.github.io/crimlab/reference/limfrom.md),
[`qltyrcd`](https://rootcoder007.github.io/crimlab/reference/qltyrcd.md),
[`dastd`](https://rootcoder007.github.io/crimlab/reference/dastd.md).

## Examples

``` r
# made-up percentages for three groups and their pooled total, with
# published limits and one estimate flagged E
d <- data.frame(
  group = c("First Nations", "Metis", "Inuit", "Total Indigenous",
            "Non-Indigenous"),
  pct = c(45.8, 52.8, 47.8, 49, 61.5),
  lo = c(43.9, 50.6, 43.1, 47.5, 61),
  hi = c(47.7, 55, 52.5, 50.5, 62),
  q = c("", "", "E", "", "")
)
indgap(d, "group", "pct", "lo", "hi", qual = "q")
#> Gap against a reference group
#> -----------------------------
#> Reference: Non-Indigenous
#> Groups:    4
#> Intervals: 95%
#> 
#>                     group  estimate  refest   ratio        ratio CI  difference   difference CI
#>             First Nations      45.8    61.5  0.7447  [0.714, 0.777]       -15.7  [-17.7, -13.7]
#>                     Metis      52.8    61.5  0.8585  [0.823, 0.896]        -8.7    [-11, -6.44]
#>                   Inuit *      47.8    61.5  0.7772  [0.704, 0.858]       -13.7  [-18.4, -8.97]
#> Total Indigenous (pooled)      49.0    61.5  0.7967  [0.772, 0.822]       -12.5  [-14.1, -10.9]
#> 
#> Below the reference (interval below 1): First Nations, Metis, Inuit * and Total Indigenous (pooled).
#> 
#> ratio: the group's estimate over the reference's; difference: the
#> group's less the reference's, in the estimate's own units.
#> (pooled): a pan-Indigenous total, shown beside, not a distinct group (poogrp).
#> Report a gap beside each group's own level, not alone.
#> * use with caution: the source flags the estimate (E, high sampling variability).
# your own table: made-up numbers in another country's layout, with relative
# standard errors instead of limits and that agency's flags
own <- data.frame(
  group = c(
    "Aboriginal",
    "Torres Strait Islander",
    "Aboriginal and Torres Strait Islander",
    "Non-Indigenous Australians"
  ),
  pct = c(31.2, 27.5, 30.8, 18.9),
  rse = c(4.1, 22.0, 3.9, 1.2),
  flag = c("", "*", "", "")
)
own <- cbind(own, limfrom(own$pct, cv = own$rse)[c("lower", "upper")])
own$q <- qltyrcd(own$flag, suppr = c("np", "**"), caution = "*")
indgap(
  own,
  "group",
  "pct",
  "lower",
  "upper",
  ref = "Non-Indigenous Australians",
  qual = "q",
  pooled = "Aboriginal and Torres Strait Islander"
)
#> Gap against a reference group
#> -----------------------------
#> Reference: Non-Indigenous Australians
#> Groups:    3
#> Intervals: 95%
#> 
#>                                          group  estimate  refest  ratio       ratio CI  difference  difference CI
#>                                     Aboriginal      31.2    18.9  1.651    [1.52, 1.8]        12.3   [9.75, 14.8]
#>                       Torres Strait Islander *      27.5    18.9  1.455  [0.917, 2.31]         8.6  [-3.27, 20.5]
#> Aboriginal and Torres Strait Islander (pooled)      30.8    18.9  1.630    [1.5, 1.77]        11.9    [9.5, 14.3]
#> 
#> Above the reference (interval above 1): Aboriginal and Aboriginal and Torres Strait Islander (pooled).
#> Not distinguishable from 1 (interval includes it): Torres Strait Islander *.
#> 
#> ratio: the group's estimate over the reference's; difference: the
#> group's less the reference's, in the estimate's own units.
#> (pooled): a pan-Indigenous total, shown beside, not a distinct group (poogrp).
#> Report a gap beside each group's own level, not alone.
#> * use with caution: the source flags the estimate (E, high sampling variability).
```
