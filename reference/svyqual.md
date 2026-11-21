# Release quality of a survey estimate

Grades an estimate by the rule of Health Canada's guides to the Canadian
surveys of alcohol and drug use (Health Canada 2022, Table 1): by the
unweighted number of observations behind it and by its coefficient of
variation. `"acceptable"`: at least `nmin` observations with the
characteristic (and, for a proportion, at least `dmin` in the domain)
and a coefficient of variation of at most `cvlo` per cent. `"marginal"`:
the same counts and a coefficient above `cvlo` up to `cvhi`; to be
released with a warning. `"unacceptable"`: fewer observations, or a
coefficient above `cvhi`; not to be released.

## Usage

``` r
svyqual(cv, n, d = NULL, cvlo = 16.6, cvhi = 33.3, nmin = 10, dmin = 50)
```

## Arguments

- cv:

  The coefficient of variation of the estimate, in per cent (as
  [`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
  reports it).

- n:

  The unweighted number of observations with the characteristic (the
  numerator of a proportion; for a mean or a total, the observations
  that contribute).

- d:

  For a proportion, the unweighted number of observations in the domain
  (the denominator); `NULL` (the default) for a mean or a total.

- cvlo, cvhi:

  The coefficients of variation, in per cent, up to which an estimate is
  acceptable and marginal.

- nmin, dmin:

  The fewest observations: with the characteristic, and in the domain.

## Value

A character vector: `"acceptable"`, `"marginal"` or `"unacceptable"`;
`NA` where the counts suffice and `cv` is missing.

## Details

The thresholds are arguments: other agencies set others, and Statistics
Canada's own have changed over the years.

Errors: `cv` not numeric or negative; `n` or `d` not numeric or
negative; a threshold that is not one number from 0, or `cvlo` above
`cvhi`; arguments of different lengths (a length of 1 is repeated).

## References

Health Canada (2022). *Canadian Postsecondary Education Alcohol and Drug
Use Survey, 2021-22: Public Use Microdata File user guide*. Section 7.3
and Table 1, quality level guidelines for estimates.

## See also

[`svymx()`](https://rootcoder007.github.io/crimlab/reference/svymx.md)
and
[`svypci()`](https://rootcoder007.github.io/crimlab/reference/svypci.md).

## Examples

``` r
svyqual(cv = c(8, 20, 40, 8, 8), n = c(120, 60, 30, 6, 40),
        d = c(900, 500, 300, 200, 45))
#> [1] "acceptable"   "marginal"     "unacceptable" "unacceptable" "unacceptable"
# a mean or a total: no denominator
svyqual(12, 35)
#> [1] "acceptable"
try(svyqual(-1, 10))
#> Error : svyqual: `cv` must be numbers from 0.
```
