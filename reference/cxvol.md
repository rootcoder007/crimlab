# Alert complexity and placement volatility, pooled or by fiscal year

The analysis of Ruhela (2026) on any release of Ontario's restrictive
confinement records: from the person-year table
([`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md))
it gives, for the years pooled or for each fiscal year, the share of
person-years with complex alert profiles (`cx = 1`), mean placement
volatility by group, the crude Poisson rate ratio, the one-to-one
nearest-neighbour propensity matching without replacement of each
complex person-year to a control on age group, gender and fiscal year
([`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md),
in 'MatchIt”s default order), the count model on the matched frame with
a random intercept per regional cluster
([`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)),
and the bounds on the average effect on the matched frame with no
assumption and under monotone treatment response and selection
([`mtrbnd()`](https://rootcoder007.github.io/crimlab/reference/mtrbnd.md)).

## Usage

``` r
cxvol(x, by = c("pool", "year"), yterm = TRUE, fam = c("nb", "poisson"), ...)
```

## Arguments

- x:

  A person-year table from
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md),
  or what
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  takes (a data frame of the release or a table key).

- by:

  `"pool"` for the years together, `"year"` for each year.

- yterm:

  `TRUE` to adjust for the fiscal year in the pooled propensity and
  count models.

- fam:

  The family of
  [`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md).

- ...:

  Passed to
  [`otpyr()`](https://rootcoder007.github.io/crimlab/reference/otpyr.md)
  when `x` is not a person-year table.

## Value

An object of class `cxvol`, a list: `res`, one row per sample (`sample`,
`n`, `n1`, `share`, `vm`, `vm1`, `vm0`, `crude`, `nm`, `vm1m`, `vm0m`,
`irr`, `lo`, `hi`, `est`, `se`, `p`, `ncl`, `fam`, `mlo`, `mhi`, `bhi`:
the Manski bounds and the upper bound under both monotonicity
assumptions, whose lower bound is 0; `ymax`, the largest count on the
matched frame, the range used by the bounds); `fits`, the
[`cntmix()`](https://rootcoder007.github.io/crimlab/reference/cntmix.md)
fits; `by`; `yterm`.

## Details

A single fiscal year cannot identify a year term: `by = "pool"` with one
year and `yterm = TRUE` is refused; `by = "year"` drops the term in each
year and says so. A sample with no complex person-year, or none without,
is refused, naming it.

## References

Ruhela, V. S. (2026). Alert complexity and placement volatility in
Ontario restrictive confinement data. Major research paper, Centre for
Criminology and Sociolegal Studies, University of Toronto.

## Examples

``` r
set.seed(4)
n <- 600
p <- data.frame(year = sample(2025:2026, n, TRUE),
                gender = sample(c("Male", "Female"), n, TRUE, c(9, 1)),
                age = sample(c("18 to 24", "25 to 49", "50+"), n, TRUE),
                rc = sample(c("A", "B", "A + B", "C"), n, TRUE),
                cx = rbinom(n, 1, 0.2))
p$vm <- rpois(n, exp(-2 + 0.3 * p$cx + (p$rc == "A + B")))
cxvol(p, fam = "poisson")
#> cxvol: alert complexity (cx) and placement volatility (vm), years pooled
#>     sample   n     cx    vm1    vm0 crude matched    irr     lo    hi     p
#>  2025-2026 600 0.1917 0.2174 0.2062 1.054     230 0.8477 0.5012 1.434 0.538
#>     bound
#>  -0.06087
#> irr: rate ratio of the matched count model (poisson); bound: upper end of [0, bound] under monotone response and
#> selection (below 0, the data contradict the two assumptions)
```
