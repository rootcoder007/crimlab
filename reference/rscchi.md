# Rao-Scott chi-square test for a survey table

The test of independence of two categorical variables in a two-way table
from a stratified, clustered and weighted sample (a victimisation or
justice survey, say), with Rao and Scott's corrections for the design,
to go with the estimates of
[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
and the tables of
[`ctabx()`](https://rootcoder007.github.io/crimlab/reference/ctabx.md).
Pearson's statistic on such data is too large: the design makes the cell
proportions vary more than simple random sampling, and the corrections
divide it by the design effects.

## Usage

``` r
rscchi(x, y, w = NULL, strata = NULL, psu = NULL, fpc = NULL)

# S3 method for class 'rscchi'
print(x, ...)
```

## Arguments

- x, y:

  The two variables: labels (numbers, text or a factor), one per
  respondent, each with two labels or more; for
  [`print()`](https://rdrr.io/r/base/print.html), `x` is an object of
  class `rscchi`.

- w:

  `NULL` or the sampling weights: positive numbers.

- strata:

  `NULL` or the stratum of each respondent: labels.

- psu:

  `NULL` (each respondent its own) or the first-stage unit of each
  respondent: labels.

- fpc:

  `NULL` or the number of first-stage units in the population of each
  respondent's stratum (constant within strata).

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `rscchi`, a list: `tab`, the weighted cell
proportions (a matrix, `x` by `y`); `test`, a data frame (`test`:
`pearson`, `rs1`, `rs2`, `f`; `stat`, `df1`, `df2`, `p`); `deff`, the
generalised design effects; `dbar` and `a2`, their mean and squared
coefficient of variation; and `des` (`n`, `strata`, `psu`, `df`).

## Details

The cell proportions `p` are weighted; Pearson's `X^2` is computed on
them scaled to the number of respondents `n`,
`X^2 = n sum (p_ij - p_i. p_.j)^2 / (p_i. p_.j)`. Their covariance `V`
comes by linearisation over the first-stage units within strata, as
[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
reads the design (Lumley 2004). The generalised design effects `deff`
are the eigenvalues of `Delta = (C' D^-1 C / n)^-1 (C' D^-1 V D^-1 C)`,
`C` a basis of the interaction contrasts and `D` the diagonal of `p`
(Rao and Scott 1984). The tests: `pearson`, `X^2` on
`K = (r - 1)(c - 1)` degrees of freedom, wrong under the design; `rs1`,
the first-order correction `X^2 / mean(deff)` on `K`; `rs2`, the
second-order (Satterthwaite) correction `X^2 / (mean(deff) (1 + a^2))`
on `K / (1 + a^2)`, `a` the coefficient of variation of the design
effects; and `f`, `X^2 / sum(deff)` as an F on `d0 = K / (1 + a^2)` and
`d0 nu` degrees of freedom, `nu` the first-stage units less the strata,
the form Thomas and Rao (1987) found to hold its level best. They agree
with `survey::svychisq()` (`"Chisq"` for `rs1`, `"F"` for `f`).

## References

Rao, J. N. K. and Scott, A. J. (1981). The analysis of categorical data
from complex sample surveys: chi-squared tests for goodness of fit and
independence in two-way tables. *Journal of the American Statistical
Association*, 76(374), 221-230.
[doi:10.1080/01621459.1981.10477633](https://doi.org/10.1080/01621459.1981.10477633)

Rao, J. N. K. and Scott, A. J. (1984). On chi-squared tests for multiway
contingency tables with cell proportions estimated from survey data.
*The Annals of Statistics*, 12(1), 46-60.
[doi:10.1214/aos/1176346391](https://doi.org/10.1214/aos/1176346391)

Thomas, D. R. and Rao, J. N. K. (1987). Small-sample comparisons of
level and power for simple goodness-of-fit statistics under cluster
sampling. *Journal of the American Statistical Association*, 82(398),
630-636.
[doi:10.1080/01621459.1987.10478476](https://doi.org/10.1080/01621459.1987.10478476)

Lumley, T. (2004). Analysis of complex survey samples. *Journal of
Statistical Software*, 9(8), 1-19.
[doi:10.18637/jss.v009.i08](https://doi.org/10.18637/jss.v009.i08)

## Examples

``` r
set.seed(1)
d <- data.frame(region = rep(c("north", "south"), each = 60),
                block = rep(1:12, each = 10),
                age = sample(c("young", "middle", "old"), 120, TRUE),
                w = runif(120, 50, 150))
d$victim <- rbinom(120, 1, ifelse(d$age == "young", 0.35, 0.2))
rscchi(d$age, d$victim, d$w, d$region, d$block)
#> rscchi: 3 x 2 table, 120 respondents, 2 strata, 12 first-stage units
#>     test  stat   df1   df2       p
#>  pearson 5.539 2.000    NA 0.06271
#>      rs1 3.254 2.000    NA 0.19650
#>      rs2 3.134 1.926    NA 0.19750
#>        f 1.627 1.926 19.26 0.22264
```
