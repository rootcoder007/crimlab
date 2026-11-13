# Oaxaca-Blinder decomposition of a gap between two groups

Splits the gap in mean outcome between two groups (a sentence length
between two groups of defendants, say) into the part their covariates'
means account for and the part they do not, by least squares in each
group (Oaxaca 1973; Blinder 1973), with delta-method standard errors
(Jann 2008).

## Usage

``` r
oaxbd(y, x, grp, ref = 1, conf = 0.95)
```

## Arguments

- y:

  The outcomes: finite numbers, one per case.

- x:

  The covariates: `NULL` (the constant only), a numeric vector, a
  numeric matrix with one row per case, or a data frame of numeric
  columns; a constant is added. Factors are refused: code them as
  indicator columns.

- grp:

  The group of each case: labels, exactly two groups, each with more
  cases than coefficients.

- ref:

  Group B, whose coefficients price the explained part of the first
  decomposition: its label, or its position from 1 among the two (in the
  order of `grp`); the other is group A.

- conf:

  The confidence level of the intervals, strictly between 0 and 1.

## Value

A list: `groups`, the labels of `a` and `b`; `means`, a data frame with
one row per term (`const` first): `term`, `a`, `b`; `coef`: `term`, `a`,
`ase`, `b`, `bse`; `dec`, one row per part (`gap`, `exp`, `unexp`,
`expa`, `unexpa`, `endow`, `coef`, `inter`): `part`, `est`, `se` and the
normal interval `lo`, `hi`; `by`: `term`, `exp`, `expse`, `unexp`,
`unexpse` and `shift`.

## Details

With `D` the difference of the covariate means (group A minus group B,
the constant first) and `Delta` the difference of the coefficients, the
gap is exactly `D'b_B + xbar_A'Delta` (explained at B's coefficients,
`exp`, and unexplained, `unexp`), and equally `D'b_A + xbar_B'Delta`
(`expa`, `unexpa`), and `D'b_B + xbar_B'Delta + D'Delta` (endowments,
coefficients and their interaction). The two explained parts differ by
exactly the interaction, so how much of a gap is "explained" depends on
which group's coefficients are taken as the standard: `oaxbd()` reports
both. `by` splits the first decomposition by covariate; its unexplained
lines are not identified one by one: moving a covariate's origin by `c`
for everyone moves `c` times its `shift` between that line and the
constant's, leaving the total alone (Oaxaca and Ransom 1999), and for a
factor coded as indicators the split depends on the category left out.
See `proof("oaxbd")`.

The standard errors treat the two groups, and the means and coefficients
within each, as independent, and leave out a term of smaller order, as
Jann (2008) does. The unexplained part is what the model's covariates do
not account for. It is not a measure of discrimination: it holds every
influence the covariates miss, and the covariates can themselves carry
earlier unequal treatment.

## References

Oaxaca, R. (1973). Male-female wage differentials in urban labor
markets. *International Economic Review*, 14(3), 693-709.
[doi:10.2307/2525981](https://doi.org/10.2307/2525981)

Blinder, A. S. (1973). Wage discrimination: reduced form and structural
estimates. *Journal of Human Resources*, 8(4), 436-455.
[doi:10.2307/144855](https://doi.org/10.2307/144855)

Oaxaca, R. L. and Ransom, M. R. (1999). Identification in detailed wage
decompositions. *Review of Economics and Statistics*, 81(1), 154-157.
[doi:10.1162/003465399767923908](https://doi.org/10.1162/003465399767923908)

Jann, B. (2008). The Blinder-Oaxaca decomposition for linear regression
models. *The Stata Journal*, 8(4), 453-479.
[doi:10.1177/1536867X0800800401](https://doi.org/10.1177/1536867X0800800401)

## Examples

``` r
set.seed(3)
g <- rep(c("a", "b"), each = 200)
sev <- rnorm(400, ifelse(g == "a", 5, 4))
pri <- rpois(400, ifelse(g == "a", 2, 1))
len <- 2 + 0.8 * sev + 0.5 * pri + ifelse(g == "a", 1, 0) + rnorm(400)
r <- oaxbd(len, data.frame(sev, pri), g, ref = "b")
r$dec
#>     part        est        se         lo        hi
#> 1    gap 2.23260343 0.1384465  1.9612533 2.5039535
#> 2    exp 1.40297787 0.1668697  1.0759193 1.7300365
#> 3  unexp 0.82962555 0.1680470  0.5002595 1.1589916
#> 4   expa 1.43218736 0.1341264  1.1693045 1.6950703
#> 5 unexpa 0.80041606 0.1344082  0.5369809 1.0638512
#> 6  endow 1.40297787 0.1668697  1.0759193 1.7300365
#> 7   coef 0.80041606 0.1344082  0.5369809 1.0638512
#> 8  inter 0.02920949 0.1641780 -0.2925734 0.3509924
```
