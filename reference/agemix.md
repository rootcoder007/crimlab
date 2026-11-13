# The age-crime curve of a mixture of offender types

The aggregate age-crime curve is the share-weighted mean of the curves
of latent types (adolescence-limited, life-course persistent, and so
on), and one aggregate curve fits many mixtures: the aggregate alone
cannot tell one type from several (Moffitt 1993; Gottfredson and Hirschi
1990). What does tell them apart is the covariance of an individual's
offending at two ages: with one type it is 0 (given the type, ages are
independent), with several it is the covariance across types of their
curves at the two ages, which `cov` gives for every pair of ages. With a
panel of individuals observed at two ages, the observed covariance and
its permutation test say whether one type is enough. See
`proof("agemix")`.

## Usage

``` r
agemix(
  curves = NULL,
  shares = NULL,
  ages = NULL,
  panel = NULL,
  nperm = 2000,
  seed = NULL
)
```

## Arguments

- curves:

  `NULL`, or the curves of the types: a numeric matrix or a data frame,
  one row per age and one column per type (at most 30).

- shares:

  The share of each type: numbers from 0 that sum to 1.

- ages:

  `NULL` (1, 2, ...) or the age of each row of `curves`.

- panel:

  `NULL`, or the outcomes of individuals at two ages: a numeric matrix
  or data frame of two columns, at least three rows.

- nperm:

  The number of permutations of the second age's outcomes for the panel
  test: a whole number from 0.

- seed:

  `NULL` (the draws follow
  [`set.seed()`](https://rdrr.io/r/base/Random.html)) or a whole number
  from 0.

## Value

An object of class `agemix`, a list: `curve` (`age`, `agg` the
aggregate, then the types' curves; with `curves`), `cov` (the covariance
across types between ages, a matrix named by age; with `curves`) and
`panel` (`cov`, `p` the permutation p-value, `NA` without permutations,
`n`, `nperm`; with `panel`); the parts not asked for are `NULL`.

## References

Moffitt, T. E. (1993). Adolescence-limited and life-course-persistent
antisocial behavior: a developmental taxonomy. *Psychological Review*,
100(4), 674-701.
[doi:10.1037/0033-295X.100.4.674](https://doi.org/10.1037/0033-295X.100.4.674)

## Examples

``` r
ages <- 12:30
al <- dnorm(ages, 17, 2)
lcp <- rep(0.08, length(ages))
r <- agemix(cbind(al = al, lcp = lcp), c(0.8, 0.2), ages)
r$cov["15", "25"]
#> [1] -0.0005241738
```
