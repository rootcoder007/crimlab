# Concentration of crime across places, people or units

How unequally a total is spread across units: the Gini coefficient and
the Lorenz curve, the share held by the top units, and what chance alone
would show. A few places holding most of the incidents, or a few people
most of the placements, is the law of crime concentration (Weisburd
2015); with sparse counts much of the Gini comes from the units with no
events at all, and a Poisson process at the same mean already leaves
many units at zero (Bernasco and Steenbeek 2017).

## Usage

``` r
cconc(
  x,
  w = NULL,
  by = NULL,
  fracs = c(0.01, 0.05, 0.1, 0.25, 0.5),
  pop = NULL,
  nsim = 2000,
  seed = NULL
)
```

## Arguments

- x:

  The count of each unit (numbers from 0 with a positive total), or the
  place of each event in the order they happened (labels: text or a
  factor).

- w:

  `NULL`, or a multiplicity for each count (positive numbers): a banded
  table of people by number of placements gives each band's count once
  with its number of people. Counts only.

- by:

  `NULL`, or a group for each unit (each event, for places): labels. The
  Gini is also given within each group.

- fracs:

  The top shares to report: shares from 0 to 1 of the units.

- pop:

  `NULL`, or the population of each unit (positive numbers; counts
  only): the urban-scaling exponent of the counts on the population is
  fitted, `log x = log pref + beta log pop`, over the units with events.

- nsim:

  The number of Poisson samples at the observed mean for the null
  distribution of the Gini (0 for none): a whole number. It needs a
  whole number of units.

- seed:

  `NULL` (the draws follow
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) or a whole number
  from 0.

## Value

An object of class `cconc`, a list: `stat`, a named vector (`units`,
`total`, `mean`, `var`, `disp`, `ivar`, `isd`, `zero`, `znull`, `zgap`,
`gini`, `gpos`, `gnull`, `vpop`, `cvb`, `nullb`); `top` (`frac`,
`units`, `share`); `lorenz` (`pop`, `share`, from the origin); `null`,
`NULL` or a named vector (`nsim`, `mean`, `q95`, `p` the share of
simulated Gini coefficients at least the observed one with the observed
counted in, and `excess` the observed minus the simulated mean); `by`,
`NULL` or a data frame (`group`, `units`, `total`, `gini`, `NA` for a
group without events); `scale`, `NULL` or a named vector (`beta`,
`pref`, `se`, `r2`, `used`, the first four `NA` with fewer than three
units with events); `grow`, `NULL` or a list (`mhat`, `slope`, `curve`
with `i`, `dist`, `expd`, `lo`, `hi`); and `places`, the labels of the
places (`NULL` for counts).

## Details

The Gini coefficient here is the population form,
`1 - sum p (L[i-1] + L[i])` over the units in ascending order (`p` a
unit's share of the units, `L` the cumulative share of the total), which
for unit weights is `2 sum i x[(i)] / (n T) - (n + 1) / n`. It splits
exactly as `zero + (1 - zero) gpos`, the share of units at zero plus the
rest times the Gini of the units with events; `gnull` is the same split
at the share of zeros a Poisson process at the observed mean leaves,
`exp(-mean)`. The variance over the mean (`disp`) is 1 for Poisson
counts, and its excess estimates the variance of the places' underlying
rates (`ivar`). The Gini of counts is at most `sqrt(vpop / 2) / mean`
(`cvb`), so for Poisson counts chance alone can produce a Gini up to
`1 / sqrt(2 mean)` (`nullb`); the simulated null gives its distribution.
See `proof("cconc")`.

Given the place of each event in order (`x` as labels), the counts are
per place, and `grow` follows the number of distinct places as events
accumulate against preferential allocation, where an event goes to a new
place with chance `M / (M + i)` (the Chinese restaurant process): `mhat`
is the `M` that matches the final count, `expd` the expected curve, `lo`
and `hi` its bounds, and `slope` the log-log slope over the second half.

## References

Weisburd, D. (2015). The law of crime concentration and the criminology
of place. *Criminology*, 53(2), 133-157.
[doi:10.1111/1745-9125.12070](https://doi.org/10.1111/1745-9125.12070)

Bernasco, W. and Steenbeek, W. (2017). More places than crimes:
implications for evaluating the law of crime concentration at place.
*Journal of Quantitative Criminology*, 33(3), 451-467.
[doi:10.1007/s10940-016-9324-7](https://doi.org/10.1007/s10940-016-9324-7)

Bettencourt, L. M. A., Lobo, J., Helbing, D., Kuhnert, C. and West, G.
B. (2007). Growth, innovation, scaling, and the pace of life in cities.
*Proceedings of the National Academy of Sciences*, 104(17), 7301-7306.
[doi:10.1073/pnas.0610172104](https://doi.org/10.1073/pnas.0610172104)

## Examples

``` r
set.seed(2)
n <- rnbinom(200, mu = 2, size = 0.4)
r <- cconc(n, nsim = 200, seed = 1)
r$stat[c("gini", "zero", "gpos", "gnull")]
#>      gini      zero      gpos     gnull 
#> 0.7151670 0.4850000 0.4469262 0.4973515 
r$top
#>   frac units      share
#> 1 0.01     2 0.08977035
#> 2 0.05    10 0.29645094
#> 3 0.10    20 0.47181628
#> 4 0.25    50 0.80584551
#> 5 0.50   100 0.99373695
# people by number of placements, from a banded table
cconc(c(1, 2, 4, 8), w = c(500, 120, 40, 9), nsim = 0)$stat["gini"]
#>      gini 
#> 0.2599236 
# the place of each incident, in order
cconc(sample(letters[1:8], 60, TRUE, prob = 8:1), nsim = 0)$grow$mhat
#> [1] 2.254821
```
