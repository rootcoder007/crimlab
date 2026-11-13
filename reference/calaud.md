# Calibration audit of a predictive-policing ranking

Does the order in which an algorithm ranks areas by predicted risk match
the order of what was later recorded there, and does the disagreement
follow the areas' demographic make-up? `calaud()` ranks the areas both
ways, measures the agreement by Spearman's correlation, and gives each
group the mean gap between its areas' two ranks, with a permutation
test. It generalises the district-level analysis of the SciencesPo
predictive-policing study of Chicago (Lacherade et al. 2021) to any
place: the five columns of the audit are `area`, `risk`, `outcome`,
`group` and `population`.

## Usage

``` r
calaud(
  area,
  risk = NULL,
  outcome = NULL,
  group = NULL,
  population = NULL,
  nperm = 999,
  seed = NULL
)
```

## Arguments

- area:

  The area of each record: labels; or a data frame with columns `area`,
  `risk`, `outcome`, `group` and optionally `population`, the other
  arguments then left out.

- risk:

  The risk the algorithm predicts for each record: finite numbers.

- outcome:

  The recorded outcome of each record (a 0 or 1, or a count): finite
  numbers.

- group:

  The group of each record (or area): labels.

- population:

  `NULL`, or the population of each area: positive, one per record (the
  same within an area) or one per area named by area.

- nperm:

  The number of permutations: a whole number from 0 (none).

- seed:

  `NULL` to draw the seed of the bootstrap (or the permutations) from
  R's generator, so that
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it, or one
  whole number from 0.

## Value

A list: `areas`, a data frame with one row per area (sorted by label):
`area`, `n` (records), `risk`, `rate`, `group`, `riskrk`, `raterk` and
`gap`; `groups`, one row per group: `group`, `n` (areas), `gap` (mean)
and `p` (`NaN` without areas or permutations); `spear`, a named vector:
`rho`, `s` (Spearman's statistic, the sum of squared rank differences),
`pex`, `t` and `pt`; `ties` (`TRUE` when either ranking has ties);
`between`, `ss` and `p`; and `nperm`.

## Details

Records are rolled up to areas first: an area's risk is the mean over
its records, its rate the mean outcome (or, with `population`, the
summed outcome per head), and its group the most frequent label among
its records (the first in order on ties). One record per area is the
per-area audit itself. Rank 1 is the highest; tied values share the mean
of their ranks. An area's gap is its outcome rank minus its risk rank:
positive when the algorithm ranks it as riskier than its recorded
outcome warrants, negative when less risky.

Spearman's rho is Pearson's correlation of the two rankings (Spearman
1904). `pex` is its two-sided p-value by algorithm AS 89 (Best and
Roberts 1975): exact for nine areas or fewer, from the Edgeworth series
up to 1289, and `NA` when either ranking has ties; `pt` is the
large-sample one from the t distribution, as
[`stats::cor.test()`](https://rdrr.io/r/stats/cor.test.html) gives them.
A group's mean gap is tested against shuffling the groups across the
areas, `nperm` times; `between` tests all groups at once by the
between-group sum of squares of the gaps (the gaps sum to 0).

The audit compares the algorithm with the records, not with crime: a
recorded outcome follows where the police looked as well as what
happened, so a group's areas can look over-predicted because their
incidents were recorded less, or under-predicted because they were
recorded more (see
[`fairsim()`](https://rootcoder007.github.io/crimlab/reference/fairsim.md)).
A gap says the two orders disagree for that group's areas; it does not
say which order is right or why.

## References

Spearman, C. (1904). The proof and measurement of association between
two things. *American Journal of Psychology*, 15(1), 72-101.
[doi:10.2307/1412159](https://doi.org/10.2307/1412159)

Best, D. J. and Roberts, D. E. (1975). Algorithm AS 89: the upper tail
probabilities of Spearman's rho. *Journal of the Royal Statistical
Society, Series C*, 24(3), 377-379.
[doi:10.2307/2347111](https://doi.org/10.2307/2347111)

Lacherade, P., Szabo, B., Krikava, J. and Aeby, M. (2021). Predictive
policing in the City of Chicago. Group project, Decoding Biases in AI,
Sciences Po. <https://github.com/marquito95/Predictive-policing-Chicago>

## Examples

``` r
set.seed(3)
d <- fairsim(4000, nareas = 24, bias = 0.5, seed = 1)
# the recorded outcome is `det`: the incidents that were detected
r <- calaud(d$area, d$risk, d$det, d$group, nperm = 499, seed = 1)
r$groups
#>   group  n gap     p
#> 1     A 12 -12 0.002
#> 2     B 12  12 0.002
r$spear
#>           rho             s           pex             t            pt 
#> -7.400000e-01  4.002000e+03  5.640811e-05 -5.160381e+00  3.575934e-05 
```
