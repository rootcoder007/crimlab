# Healing as a process: change in published indicators between two periods

For each group and indicator, the change between two survey periods with
its confidence interval, wellness indicators first: change over time
read as Indigenous scholars describe mental health, as a process of
healing, not as a fixed deficit.

## Usage

``` r
healchg(
  data,
  grp,
  ind,
  prd,
  est,
  lower,
  upper,
  from,
  to,
  well = NULL,
  qual = NULL,
  conflvl = 0.95,
  distress = NULL
)
```

## Arguments

- data:

  A data frame, one row per group, indicator and period.

- grp, ind, prd, est, lower, upper:

  Column names: the change is in the estimate's own units (percentage
  points for percentages).

- from, to:

  The two periods.

- well, distress:

  Optional indicators of wellness (listed first) and of distress; with
  `distress` given, any other indicator is marked `"other"`.

- qual:

  Optional column of quality flags: `"E"` rows are marked use with
  caution, and suppressed ones (`"F"`, `"x"`, `".."`) get no estimate.

- conflvl:

  Confidence level of each change's interval (default 0.95).

## Value

A data frame of class `healchg`: `group`, `indicator`, `side`, `from`,
`to`, `change`, `lower`, `upper`, `usewithcau`.

## Details

Stewart (2008) describes Indigenous mental health as wellness and as a
process of healing; a process is seen in change. The change \\p_2 -
p_1\\ has the interval \\p_2 - p_1 \pm z\sqrt{s_1^2 + s_2^2}\\, with the
standard errors from the published intervals and the two periods'
samples independent (separate cycles of a survey are). A change is only
as comparable as the two periods' questions and designs: Statistics
Canada flags breaks in a series, and a change across one is not a change
in the people.

## Acknowledgement

This estimator puts to quantitative use the framework of Professor
Suzanne Stewart (Yellowknives Dene First Nation, University of Toronto).

See the acknowledgements in the section on Indigenous health of the
package's vignette
([`vignette("crimlab")`](https://rootcoder007.github.io/crimlab/articles/crimlab.md)).

## References

Stewart SL (2008), see
[`holbal`](https://rootcoder007.github.io/crimlab/reference/holbal.md).

## See also

[`welprof`](https://rootcoder007.github.io/crimlab/reference/welprof.md).

## Examples

``` r
# your own table: made-up percentages for one community in two surveys
own <- data.frame(
  community = "Lakeside",
  item = rep(c("Good mental health", "Anxiety"), 2),
  survey = rep(c("2018", "2023"), each = 2),
  pct = c(55, 14, 61, 11),
  lo = c(50, 11, 56, 8),
  hi = c(60, 17, 66, 14)
)
healchg(own, "community", "item", "survey", "pct", "lo", "hi",
        from = "2018", to = "2023", well = "Good mental health",
        distress = "Anxiety")
#> Healing as change between two periods
#> -------------------------------------
#> Periods:    2018 to 2023
#> Groups:     Lakeside
#> Indicators: 2
#> Units:      the estimate's own (percentage points for percentages)
#> 
#> Good mental health (wellness)
#>    group  from  to  change       interval
#> Lakeside    55  61       6  [-1.07, 13.1]
#> 
#> Anxiety (distress)
#>    group  from  to  change       interval
#> Lakeside    14  11      -3  [-7.24, 1.24]
#> 
#> 0 of 2 changes have an interval clear of zero.
#> Wellness: 0 rose and 0 fell by more than chance; distress: 0 rose and 0 fell.
#> 
#> Wellness indicators come first: a change is read as part of healing.
#> A change is only as comparable as the two periods' questions and
#> designs: across a break in the series it is not a change in people.
```
