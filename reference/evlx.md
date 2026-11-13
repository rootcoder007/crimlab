# E-values: how strong an unmeasured confounder must be

The E-value of VanderWeele and Ding (2017): the smallest strength of
association, on the risk-ratio scale, that an unmeasured confounder
would need with both the treatment and the outcome to explain away an
estimate; and that of the confidence limit nearer the null.

## Usage

``` r
evlx(est, lo = NULL, hi = NULL, type = "rr", rare = FALSE, sd = NULL)
```

## Arguments

- est:

  The estimate: one positive ratio, or one standardised difference.

- lo, hi:

  `NULL` or the confidence limits of a ratio.

- type:

  `"rr"`, `"or"`, `"hr"` or `"smd"`.

- rare:

  `TRUE` when the outcome is rare (an odds or hazard ratio is then read
  as a risk ratio).

- sd:

  `NULL` or the standard error of a standardised difference.

## Value

A named numeric vector: `rr` (the estimate as a risk ratio), `evalue`
and `evlim` (the limit's E-value; `NA` without limits).

## Details

For a risk ratio `RR >= 1` the E-value is `RR + sqrt(RR (RR - 1))`; a
ratio below 1 is inverted first. An odds or hazard ratio is a risk ratio
when the outcome is rare (`rare = TRUE`); otherwise an odds ratio is
taken to `sqrt(OR)` and a hazard ratio to
`(1 - 0.5^sqrt(HR)) / (1 - 0.5^sqrt(1 / HR))`. A standardised mean
difference `d` is taken to `exp(0.91 d)`, its limits to
`exp(0.91 d -+ 1.78 sd)` from its standard error `sd`. The limit's
E-value is 1 when the interval holds the null. These are the conversions
of the 'EValue' package.

## References

VanderWeele, T. J. and Ding, P. (2017). Sensitivity analysis in
observational research: introducing the E-value. *Annals of Internal
Medicine*, 167(4), 268-274.
[doi:10.7326/M16-2607](https://doi.org/10.7326/M16-2607)

## Examples

``` r
evlx(3.9, 1.8, 8.7)
#>       rr   evalue    evlim 
#> 3.900000 7.263034 3.000000 
evlx(0.8, 0.71, 0.91, type = "hr")
#>        rr    evalue     evlim 
#> 0.8567726 1.6088913 1.3360826 
evlx(0.5, type = "smd", sd = 0.1)
#>       rr   evalue    evlim 
#> 1.576173 2.529142 1.968037 
```
