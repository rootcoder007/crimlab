# Concentration index of exposure by income

Measures whether exposure to pollution is concentrated among
lower-income or higher-income people. The index runs from -1 (all the
exposure borne by the poorest) through 0 (exposure shared evenly across
income) to +1 (all of it borne by the richest); environmental-justice
analyses report it next to the burden.

## Usage

``` r
expconcidx(data, expo, income)

# S3 method for class 'expconcidx'
print(x, ...)
```

## Arguments

- data:

  A data frame with one row per person or area.

- expo, income:

  Column names (single strings): the exposure (a non-negative
  concentration) and the income measure (numeric).

- x:

  A result.

- ...:

  Not used.

## Value

A list of class `expconcidx` with `concidx` (a number in \\\[-1, 1\]\\),
`interpretation` (a sentence), `expomean`, `citation` and `extra`, a
list with `n` (rows used), `cov_h_R` (the covariance) and `curve`, a
data frame of the concentration curve: `population`, the cumulative
population share, poorest first, and `exposure`, the cumulative share of
total exposure.

## Details

Wagstaff, Paci and van Doorslaer (1991): \$\$CI = (2 / \mu) \\
\mathrm{cov}(h_i, R_i)\$\$ with \\R_i = (rank_i - 0.5) / n\\ the
fractional income rank (tied incomes share their average rank, so the
index does not depend on row order) and the population covariance.
Negative values mean lower-income units bear more exposure. Both columns
must be numeric: a factor or character income is refused rather than
ranked by its level codes, and exposures must be non-negative.

Rows with a missing exposure or income are dropped. The interpretation
calls an index below -0.01 a burden on the poor, above 0.01 a burden on
the rich, and anything between an even distribution. Income only needs
an order: an amount, a quintile number or a rank all give the same
index.

## References

Wagstaff, A., Paci, P. and van Doorslaer, E. (1991). On the measurement
of inequalities in health. Social Science and Medicine 33(5), 545-557.

## Examples

``` r
d <- data.frame(exposure = c(30, 25, 20, 15, 10), income = 1:5)
ci <- expconcidx(d, "exposure", "income")
ci$concidx
#> [1] -0.2
ci$interpretation
#> [1] "CI = -0.2000. Pro-poor exposure burden: lower-income individuals bear disproportionately higher pollution."
try(expconcidx(data.frame(exposure = 1:3, income = c("a", "b", "c")),
               "exposure", "income"))
#> Error : expconcidx: `income` must be numeric, got character.
```
