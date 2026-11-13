# Change point: a step in a series and Pettitt's test

Asks whether a series moved to a new level, and where. The step scan
takes every split after period `k` with at least `minseg` periods on
each side, and the largest between-segment sum of squares \\T_k = k (n -
k) / n (\bar y_a - \bar y_b)^2\\; its p-value comes from the same
maximum over reorderings of the series (every ordering when `n!` is at
most `nperm`, giving an exact p-value; otherwise `nperm` random ones,
`(1 + #{T* >= T}) / (1 + nperm)`), so it is a test of a step anywhere,
not of the best split alone. Pettitt's rank test takes the split with
the largest \\\|U_t\|\\, \\U_t = \sum\_{i \le t \< j} sign(y_j - y_i)\\,
with its approximate p-value \\2 \exp(-6 K^2 / (n^3 + n^2))\\. Both
assume exchangeable periods under no change: a steady trend can read as
a step.

## Usage

``` r
stpcp(y, x = NULL, minseg = 2, nperm = 9999, seed = NULL)
```

## Arguments

- y:

  The series, in period order: finite numbers.

- x:

  `NULL`, or the periods (any vector, one per value): `at` then names
  the last period before each change.

- minseg:

  The fewest periods on either side of a step: one whole number from 1.

- nperm:

  The reorderings of the step scan's null: one whole number from 1.

- seed:

  `NULL` to draw the seed from R's generator, or one whole number from
  0.

## Value

An object of class `stpcp`, a list: `tab`, one row for each test
(`test`, `k`, the last period before the change, the means `before` and
`after`, the statistic `stat` (T or K) and `p`, with `at` from `x`);
`nperm`, the orderings used; `exact`; and `n`.

## References

Pettitt, A. N. (1979). A non-parametric approach to the change-point
problem. *Journal of the Royal Statistical Society, Series C*, 28(2),
126-135. [doi:10.2307/2346729](https://doi.org/10.2307/2346729)

## Examples

``` r
stpcp(c(100, 104, 98, 60, 63, 58), x = 2018:2023)
#> stpcp: 6 periods; step p-value exact over 720 orderings
#>     test k   before    after     stat         p   at
#>     step 3 100.6667 60.33333 2440.167 0.1000000 2020
#>  pettitt 3 100.6667 60.33333    9.000 0.2907114 2020
stpcp(c(100, 104, 98, 101, 97, 60, 63, 58, 61, 59), nperm = 999, seed = 1)
#> stpcp: 10 periods; step p-value from 999 orderings
#>     test k before after   stat         p
#>     step 5    100  60.2 3960.1 0.0100000
#>  pettitt 5    100  60.2   25.0 0.0661425
```
