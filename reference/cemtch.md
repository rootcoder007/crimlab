# Coarsened exact matching

The matching of Iacus, King and Porro (2012): each covariate is
coarsened into bins, the rows are grouped by their bins on every
covariate at once, and only the groups (strata) holding both treated and
control rows are kept, so that a later regression compares like with
like without extrapolating.

## Usage

``` r
cemtch(data, d, x, cut = NULL, lcut = NULL)

# S3 method for class 'cemtch'
print(x, ...)
```

## Arguments

- data:

  A data frame.

- d:

  The name of the treatment column: 0 or 1 (or `FALSE`, `TRUE`).

- x:

  The names of the covariate columns to match on.

- cut:

  `NULL` or a named list: for a numeric covariate, the number of bins (a
  whole number from 1) or the breaks (increasing, covering its range).

- lcut:

  The same for the coarsening on which `l1` is measured.

- ...:

  Not used.

## Value

An object of class `cemtch`, a list: `keep` (logical, the rows matched),
`strata` (the stratum of each row, from 1; `NA` when not matched), `w`
(the weights, 0 when not matched), `tab` (rows by group, `all`,
`matched` and `unmatched`), `l1` (`before`, `after`), `nstrata` (the
strata kept and in all) and `cut` (the breaks used).

## Details

A numeric covariate is cut into equal-width bins over its range,
[`nclass.Sturges()`](https://rdrr.io/r/grDevices/nclass.html) of them
unless `cut` gives the number of bins or the breaks (bins closed on the
left, the last closed on both sides); a character, factor or logical
covariate is matched exactly on its values. A kept stratum `s` with
`m_T^s` treated and `m_C^s` control rows gives its treated rows weight 1
and its controls weight `(m_T^s / m_C^s) (m_C / m_T)`, `m_T` and `m_C`
the matched totals: with these weights the controls stand for the
treated rows' distribution, and a weighted regression of the outcome on
the treatment (with `w` as weights, on the rows of `keep`) estimates the
effect on the treated (Iacus, King and Porro 2012, sec. 3; Blackwell et
al. 2009, p. 529-531).

`l1` is the multivariate imbalance `L1` of Iacus, King and Porro (2011):
half the sum over the cells of a coarsening of the absolute differences
between the treated and the control shares of each cell, 0 for identical
histograms and 1 for disjoint ones, before matching and after (with the
weights). On the matching coarsening it is 0 after matching by
construction, so it is measured on the coarsening of `lcut`, by default
Scott's rule
([`nclass.scott()`](https://rdrr.io/r/grDevices/nclass.html) equal-width
bins) for each numeric covariate.

Sampson and Winter (2018) define children with a mean blood lead level
of 6 micrograms per decilitre or more as treated, match them on
categorical traits of the child, the caregiver and the neighbourhood
(Blackwell et al. 2009), report `L1` falling "from .71 to nearly 0" (p.
283), and then fit ordinary least squares on the matched sample; they
assess mediators (impulsivity, anxiety) by adding them to that
regression (p. 289).

Talayero et al. (2023) review seventeen studies of measured lead (blood,
bone or dentine) and later delinquent or criminal behaviour, rating each
for confounding, and find "limited data at the individual level"
(abstract).

## References

Iacus, S. M., King, G. and Porro, G. (2012). Causal inference without
balance checking: coarsened exact matching. *Political Analysis*, 20,
1-24. [doi:10.1093/pan/mpr013](https://doi.org/10.1093/pan/mpr013)

Iacus, S. M., King, G. and Porro, G. (2011). Multivariate matching
methods that are monotonic imbalance bounding. *Journal of the American
Statistical Association*, 106, 345-361.
[doi:10.1198/jasa.2011.tm09599](https://doi.org/10.1198/jasa.2011.tm09599)

Blackwell, M., Iacus, S., King, G. and Porro, G. (2009). cem: coarsened
exact matching in Stata. *The Stata Journal*, 9, 524-546.

Sampson, R. J. and Winter, A. S. (2018). Poisoned development: assessing
childhood lead exposure as a cause of crime in a birth cohort followed
through adolescence. *Criminology*, 56(2), 269-301.
[doi:10.1111/1745-9125.12171](https://doi.org/10.1111/1745-9125.12171)

Talayero, M. J., Robbins, C. R., Smith, E. R. and Santos-Burgoa, C.
(2023). The association between lead exposure and crime: a systematic
review. *PLOS Global Public Health*, 3(8), e0002177.
[doi:10.1371/journal.pgph.0002177](https://doi.org/10.1371/journal.pgph.0002177)

## Examples

``` r
set.seed(2)
n <- 300
d <- data.frame(pov = rbinom(n, 1, 0.4), age = runif(n, 14, 18))
d$lead <- rbinom(n, 1, plogis(-1 + 1.5 * d$pov))
d$y <- 0.5 * d$lead + d$pov + 0.2 * d$age + rnorm(n)
m <- cemtch(d, "lead", c("pov", "age"), cut = list(age = 4))
m
#> cemtch: coarsened exact matching, 8 of 8 strata kept
#>           control treated
#> all           177     123
#> matched       177     123
#> unmatched       0       0
#> L1 imbalance: 0.416 before, 0.142 after
stats::lm(y ~ lead, d, weights = m$w, subset = m$keep)
#> 
#> Call:
#> stats::lm(formula = y ~ lead, data = d, subset = m$keep, weights = m$w)
#> 
#> Coefficients:
#> (Intercept)         lead  
#>      3.9217       0.5045  
#> 
```
