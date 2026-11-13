# Inference for a difference in differences with few clusters

Tests of the two-way fixed-effects effect of
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
that stay accurate when a policy was assigned to a handful of clusters
(districts, states), where the cluster-robust t test over-rejects: the
restricted wild cluster bootstrap, and the permutation of the adoption
dates across clusters; and the power of a design simulated on real
untreated data.

## Usage

``` r
wcbx(
  y,
  unit,
  time,
  trt = NULL,
  adj = NULL,
  clus = NULL,
  meth = "wild",
  reps = NULL,
  wtype = "rad",
  ties = "count",
  stat = "t",
  ntrt = 1,
  adopt = NULL,
  eff = NULL,
  alpha = 0.05,
  power = 0.8,
  seed = NULL,
  level = 0.95
)
```

## Arguments

- y:

  The outcome: numbers, one per row; for `meth = "perm"` also a matrix
  or data frame of several outcomes, one column each.

- unit:

  The unit of each row: labels.

- time:

  The period of each row: numbers (a year, a month index); each unit at
  most once per period.

- trt:

  The treatment of each row, 0 or 1 (`"wild"`, `"perm"`; for `"perm"` it
  must stay on, from one period within each cluster); `NULL` for
  `"power"`.

- adj:

  `NULL` or time-varying controls: a data frame or list of columns, a
  numeric matrix, or one vector, one value per row. Factor, text and
  logical columns enter as indicators of each level but the first.

- clus:

  `NULL` (cluster by unit) or the cluster of each row.

- meth:

  `"wild"`, `"perm"` or `"power"`.

- reps:

  The draws: a whole number from 9; `NULL` for 9999 (200 for `"power"`).

- wtype:

  `"rad"` (Rademacher) or `"webb"`.

- ties:

  `"count"` or `"strict"`.

- stat:

  `"t"` or `"est"`: the statistic permuted.

- ntrt:

  The clusters treated in each power draw.

- adopt:

  The adoption periods drawn from, periods of `time`.

- eff:

  `NULL` (0 to two standard deviations of `y` in five steps) or the
  effects tried.

- alpha:

  The level of the power test.

- power:

  The power of the minimum detectable effect.

- seed:

  The seed of the draws: `NULL` to draw it from R's generator, or a
  whole number from 0.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

An object of class `wcbx`, a list: `res`, for `"wild"` one row (`est`,
`stat`, `p`, `lo`, `hi`), for `"perm"` one row per outcome (`outcome`,
`est`, `stat`, `p`, `pstep`), for `"power"` one row per effect (`eff`,
`power`); `meth`; `ncl`; `draws`; `enum` (`TRUE` when every draw was
enumerated); `mde` and `size` (`"power"`, else `NA`); and `level`.

## Details

`meth = "wild"` is the restricted wild cluster bootstrap-t (WCR) of
Cameron, Gelbach and Miller (2008): the model is fitted with the null
`b = b0` imposed, each bootstrap outcome is that fit plus its residuals
multiplied cluster by cluster by a random sign (`wtype = "rad"`) or one
of Webb's six values (`"webb"`), and the CR1 t statistic of
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
is recomputed. The p-value is the share of draws with `|t*|` at least
`|t|` (`ties = "count"`: the sign patterns that reproduce the data are
among the draws, as the observed assignment is in a permutation test) or
above it (`"strict"`, as 'fwildclusterboot' counts). With Rademacher
weights and `2^G <= reps` every sign pattern is enumerated. The interval
collects the `b0` not rejected at `1 - level`, by bisection with the
same draws; with Rademacher weights the smallest attainable p-value is
`2 / 2^G`, so with four clusters nothing is rejected at 5% and the
interval is the whole line: the honest answer.

`meth = "perm"` shuffles the adoption period of each cluster (all its
units adopt together; never is a period too) among the clusters, keeping
the number adopting at each date, and recomputes the t statistic
(`stat = "t"`) or the estimate (`"est"`). When the distinct assignments
number at most `reps` they are all enumerated and `p` is the exact share
at least as extreme; otherwise `p = (b + 1) / (R + 1)` over `R` random
ones. The test is Fisher's: valid in any sample when the dates were
assigned at random. With several outcomes (the columns of `y`), `pstep`
holds the Westfall-Young stepdown p-values of
[`wystep()`](https://rootcoder007.github.io/crimlab/reference/wystep.md)
over the same assignments.

`meth = "power"` draws `ntrt` clusters, gives each an adoption period
drawn from `adopt`, adds each effect of `eff` to their outcomes from
adoption on, and counts the CR1 t tests that reject at `alpha`; on real
untreated data the serial correlation and clustering are carried, as
formula-based power calculations do not (Burlig, Preonas and Woerman
2020). `size` is the rejection rate at effect 0 (it should be near
`alpha`); `mde`, where the power curve first crosses `power`, by linear
interpolation.

## References

Cameron, A. C., Gelbach, J. B. and Miller, D. L. (2008). Bootstrap-based
improvements for inference with clustered errors. *Review of Economics
and Statistics*, 90(3), 414-427.
[doi:10.1162/rest.90.3.414](https://doi.org/10.1162/rest.90.3.414)

Webb, M. D. (2023). Reworking wild bootstrap-based inference for
clustered errors. *Canadian Journal of Economics*, 56(3), 839-858.
[doi:10.1111/caje.12661](https://doi.org/10.1111/caje.12661)

Burlig, F., Preonas, L. and Woerman, M. (2020). Panel data and
experimental design. *Journal of Development Economics*, 144, 102458.
[doi:10.1016/j.jdeveco.2020.102458](https://doi.org/10.1016/j.jdeveco.2020.102458)

## Examples

``` r
d <- expand.grid(year = 2015:2024, town = LETTERS[1:8])
i <- match(d$town, LETTERS)
d$bylaw <- as.integer(i <= 3 & d$year >= 2019)
d$rate <- 30 + i + 0.5 * (d$year - 2015) - 2 * d$bylaw + sin(seq_len(80))
# eight towns: every one of the 2^8 sign patterns
wcbx(d$rate, d$town, d$year, d$bylaw)
#> wcbx: wild, 8 clusters, 256 draws (all enumerated)
#>        est      stat         p        lo        hi
#>  -1.987082 -9.207341 0.0078125 -2.524054 -1.431829
# every assignment of the by-law to three of the eight towns
wcbx(d$rate, d$town, d$year, d$bylaw, meth = "perm")$res
#>   outcome       est      stat          p      pstep
#> 1      y1 -1.987082 -9.207341 0.01785714 0.01785714
```
