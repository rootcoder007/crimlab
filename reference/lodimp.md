# Values below the limit of detection

Substitutes the values of an assay that fall below the limit of
detection (left-censored) by a fraction of that limit: `lod / 2`,
`lod / sqrt(2)` (Hornung and Reed 1990) or `lod` itself. A value not
measured (`NA`) stays `NA` and is not counted as censored. Substitution
is the simple choice; with much censoring it biases the mean and the
spread.

## Usage

``` r
lodimp(x, lod, meth = "half")
```

## Arguments

- x:

  The measured values: numbers, `NA` for not measured. A `clpre` list
  from
  [`clprep()`](https://rootcoder007.github.io/crimlab/reference/clprep.md)
  is taken too (the `tox-lodimp` preset): its elements fill the
  arguments by name.

- lod:

  The limit of detection: one positive number, or one per value (each
  sample's own limit).

- meth:

  `"half"`, `"sqrt2"` or `"lod"`.

## Value

A list: `tab`, a data frame with one row per value (`x`, `lod`, `imp`
the value after substitution, `cens` whether it was below the limit);
`frac`, the share of the measured values below the limit (`NA` when none
was measured); `nmiss`, the values not measured.

## References

Hornung, R. W. and Reed, L. D. (1990). Estimation of average
concentration in the presence of nondetectable values. *Applied
Occupational and Environmental Hygiene*, 5(1), 46-51.
[doi:10.1080/1047322X.1990.10389587](https://doi.org/10.1080/1047322X.1990.10389587)

## Examples

``` r
lodimp(c(0.4, 0.05, NA, 1.2, 0.08), lod = 0.1)
#> $tab
#>      x lod  imp  cens
#> 1 0.40 0.1 0.40 FALSE
#> 2 0.05 0.1 0.05  TRUE
#> 3   NA 0.1   NA FALSE
#> 4 1.20 0.1 1.20 FALSE
#> 5 0.08 0.1 0.05  TRUE
#> 
#> $frac
#> [1] 0.5
#> 
#> $nmiss
#> [1] 1
#> 
```
