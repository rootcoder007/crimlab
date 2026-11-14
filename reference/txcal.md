# Calibration line of an assay

Fits the calibration line of a quantitative assay (the response, an
analyte to internal standard peak-area ratio say, against the known
concentration of each calibrator) by weighted least squares, and reads
unknown samples off it. Weights `1/x^2` or `1/x` keep the precise low
calibrators from being swamped by the high ones, as is usual when the
error grows with the concentration. The limits of detection and of
quantitation follow the calibration-curve rule of ICH Q2: `k[1]` and
`k[2]` (3.3 and 10) times the residual standard error over the slope.

## Usage

``` r
txcal(conc, resp, w = "1/x^2", unk = NULL, k = c(3.3, 10))
```

## Arguments

- conc:

  The concentrations of the calibrators: positive numbers, at least
  three.

- resp:

  The response of each calibrator: finite numbers.

- w:

  The weights: `"1/x^2"`, `"1/x"`, `"none"`, or positive numbers, one
  per calibrator.

- unk:

  `NULL`, or the responses of unknown samples: numbers or `NA`.

- k:

  The multipliers of the limits of detection and quantitation: two
  positive numbers.

## Value

A list: `fit`, a named vector (`int` and `slope` of the line, `r2` the
weighted R squared, `sigma` the residual standard error on `n - 2`
degrees of freedom, `lod` and `loq` the limits); `unk`, `NULL` or a data
frame with one row per unknown (`resp`, `conc` the concentration read
off the line, `flag`: `below_lod`, `below_loq` or `quantifiable`); `w`,
the weighting.

## References

ICH (2023). *ICH Harmonised Guideline Q2(R2): Validation of Analytical
Procedures*, final version adopted 1 November 2023, section 3.2.3.3.
International Council for Harmonisation of Technical Requirements for
Pharmaceuticals for Human Use.
<https://www.ich.org/page/quality-guidelines>

## Examples

``` r
conc <- c(5, 10, 25, 50, 100, 250)
resp <- c(0.051, 0.098, 0.26, 0.49, 1.02, 2.47)
txcal(conc, resp, unk = c(0.03, 0.6))
#> $fit
#>          int        slope           r2        sigma          lod          loq 
#> 0.0005578058 0.0100118968 0.9990611257 0.0002794509 0.0921092214 0.2791188527 
#> 
#> $unk
#>   resp      conc         flag
#> 1 0.03  2.940721 quantifiable
#> 2 0.60 59.872990 quantifiable
#> 
#> $w
#> [1] "1/x^2"
#> 
```
