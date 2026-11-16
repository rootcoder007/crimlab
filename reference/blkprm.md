# Permutation test within blocks

The randomisation test of a treatment effect in a blocked (or
stratified, or matched) experiment: the treatment labels are permuted
within each block, so that every arrangement keeps each block's number
of treated units, as the design randomised them, and the observed
statistic is compared with its values over the arrangements (Rosenbaum
2002, ch. 2; Good 2005, sec. 3.6).

## Usage

``` r
blkprm(
  data,
  y,
  d,
  block,
  stat = c("diff", "strat"),
  alt = c("two.sided", "greater", "less"),
  nperm = 2000,
  exact = FALSE,
  seed = NULL
)

# S3 method for class 'blkprm'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  The name of the outcome column (numbers).

- d:

  The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`), both
  present.

- block:

  The name of the block column (labels).

- stat:

  `"diff"` (the difference in means) or `"strat"` (the block-weighted
  mean of the blocks' differences).

- alt:

  The alternative: `"two.sided"`, `"greater"` (the treated higher) or
  `"less"`.

- nperm:

  The number of arrangements drawn: a whole number from 1 to 10000000
  (not used with `exact`).

- exact:

  `TRUE` to enumerate every arrangement.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647; not used with `exact`.

- x:

  An object of class `blkprm`.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `blkprm`, a list: `est` (the observed statistic),
`p`, `nperm` (the arrangements used), `nge` (those at least as extreme),
`nblock`, `stat`, `alt`, `exact`, `dist` (the permuted values) and
`seed` (`NA` with `exact`).

## Details

`stat = "diff"` is the difference between the mean outcome of the
treated and of the controls over all units; `stat = "strat"` is the mean
over blocks of the blocks' own differences, weighted by the blocks'
sizes, the stratified estimate of the average effect (every block then
needs treated and control units). With `exact = TRUE` every arrangement
is enumerated (at most 1000000) and `p` is the share at least as extreme
as the observed value; otherwise `nperm` arrangements are drawn from the
package's generator and `p = (1 + k) / (1 + nperm)`, `k` of them at
least as extreme, which is never 0 (Phipson and Smyth 2010). Values
within `1e-9 (|T| + max|y|)` of the observed `T` count as at least as
extreme, so that arrangements equal to it in exact arithmetic are never
lost to rounding.

## References

Rosenbaum, P. R. (2002). *Observational Studies*, 2nd ed. Springer.
[doi:10.1007/978-1-4757-3692-2](https://doi.org/10.1007/978-1-4757-3692-2)

Good, P. I. (2005). *Permutation, Parametric, and Bootstrap Tests of
Hypotheses*, 3rd ed. Springer.
[doi:10.1007/b138696](https://doi.org/10.1007/b138696)

Phipson, B. and Smyth, G. K. (2010). Permutation p-values should never
be zero: calculating exact p-values when permutations are randomly
drawn. *Statistical Applications in Genetics and Molecular Biology*,
9(1), 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

## Examples

``` r
set.seed(4)
d <- data.frame(site = rep(1:6, each = 4), t = rep(c(1, 0, 1, 0), 6))
d$y <- d$site + 0.8 * d$t + rnorm(24, sd = 0.5)
blkprm(d, "y", "t", "site", nperm = 999, seed = 1)
#> blkprm: difference in means 0.9952, p = 0.001 (two.sided, 999 arrangements within 6 blocks)
blkprm(d[1:12, ], "y", "t", "site", stat = "strat", exact = TRUE)$p
#> [1] 0.01851852
```
