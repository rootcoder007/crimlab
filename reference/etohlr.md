# Ethanol in a decomposed body: drunk before death or made after

Microbes make ethanol in a body after death, and with it the higher
alcohols n-propanol and n-butanol; ethyl glucuronide, a metabolite, is
made only by a living body. So for each case: ethyl glucuronide present
(above 0) says the ethanol was taken before death; otherwise a higher
alcohol above `thr` says it was (at least partly) made after death;
otherwise the signals do not decide.

## Usage

``` r
etohlr(eth, prop = 0, but = 0, etg = NA, thr = 0)
```

## Arguments

- eth:

  The ethanol concentration: numbers from 0, one per case.

- prop, but:

  The n-propanol and n-butanol concentrations: numbers or `NA` (not
  measured), one per case or one for all.

- etg:

  The ethyl glucuronide concentration: numbers or `NA` (not measured),
  one per case or one for all.

- thr:

  The level a higher alcohol must exceed: one number from 0.

## Value

A data frame with one row per case: `eth`, `prop`, `but`, `etg`,
`higher` (the larger higher alcohol, `NA` when neither was measured) and
`verdict`, a factor (`antemortem`, `postmortem`, `indeterminate`).

## References

Kugelberg, F. C. and Jones, A. W. (2007). Interpreting results of
ethanol analysis in postmortem specimens: a review of the literature.
*Forensic Science International*, 165(1), 10-29.
[doi:10.1016/j.forsciint.2006.05.004](https://doi.org/10.1016/j.forsciint.2006.05.004)

## Examples

``` r
etohlr(c(0.8, 0.5, 0.3), prop = c(0, 0.02, NA), etg = c(1.1, NA, NA))
#>   eth prop but etg higher       verdict
#> 1 0.8 0.00   0 1.1   0.00    antemortem
#> 2 0.5 0.02   0  NA   0.02    postmortem
#> 3 0.3   NA   0  NA   0.00 indeterminate
```
