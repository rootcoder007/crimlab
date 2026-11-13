# Covariate balance and overlap

The standardised mean difference and variance ratio of each covariate
between the treated and the controls, before and after weighting, with a
verdict by a threshold; and the overlap of the propensity between the
arms.

## Usage

``` r
smdbal(x, d, w = NULL, ps = NULL, thr = 0.1)
```

## Arguments

- x:

  The covariates: a numeric matrix, a data frame or list of numeric
  columns, or one numeric vector, one row per unit.

- d:

  The treatment: 0 or 1 (or `FALSE` or `TRUE`), two units in each arm at
  least.

- w:

  `NULL` or weights from 0, one per unit (from
  [`cbatt()`](https://rootcoder007.github.io/crimlab/reference/cbatt.md),
  say).

- ps:

  `NULL` or the propensity of each unit, strictly between 0 and 1.

- thr:

  The threshold of the verdict: one positive number.

## Value

An object of class `smdbal`, a list: `tab`, one row per covariate
(`var`, `m1`, `m0`, `smd`, `vr`, and with `w` also `m1w`, `m0w`, `smdw`,
`vrw`); `overlap`, one row (`min1`, `max1`, `min0`, `max0`, `out`,
`wmax`); `bal`; `ovl`; and `thr`.

## Details

The standardised mean difference is
`(m1 - m0) / sqrt((s1^2 + s0^2) / 2)` (Imbens and Rubin 2015): means
weighted by `w` when given, the standard deviations always those of the
unweighted arms, as 'cobalt' computes it by default, so that weighting
moves only the numerator. The variance ratio is `s1^2 / s0^2`, weighted
when `w` is given (with the unbiased correction for weights). `bal` is
`TRUE` when every difference (weighted, with `w`) is below `thr` in
absolute value; 0.1 is the usual rule.

Overlap (Cole and Hernan 2008): the propensity's range in each arm, the
share of units outside the range both arms share (`out`), and the
largest inverse probability weight (`wmax`); `ovl` is `TRUE` when no
unit lies outside the common range. The propensity is the logistic
regression of `d` on `x` unless `ps` gives it.

## References

Imbens, G. W. and Rubin, D. B. (2015). *Causal Inference for Statistics,
Social, and Biomedical Sciences*. Cambridge University Press.
[doi:10.1017/CBO9781139025751](https://doi.org/10.1017/CBO9781139025751)

Cole, S. R. and Hernan, M. A. (2008). Constructing inverse probability
weights for marginal structural models. *American Journal of
Epidemiology*, 168(6), 656-664.
[doi:10.1093/aje/kwn164](https://doi.org/10.1093/aje/kwn164)

## Examples

``` r
set.seed(1)
x <- data.frame(age = rnorm(200, 40, 10), prior = rpois(200, 2))
d <- rbinom(200, 1, plogis((x$age - 40) / 10))
smdbal(x, d)
#> smdbal: balanced at 0.1: FALSE; overlap: FALSE
#>    var        m1        m0         smd        vr
#>    age 44.748673 36.687615  0.95167893 1.4812747
#>  prior  2.010989  2.119266 -0.07443683 0.7480863
#>        min1      max1       min0      max0  out     wmax
#>  0.07702527 0.9358613 0.04919905 0.7742637 0.09 12.98275
w <- cbatt(rnorm(200), d, x, est = "ebal", tgt = "att", reps = 0)$w
smdbal(x, d, w = w)$tab
#>     var        m1        m0         smd        vr       m1w       m0w
#> 1   age 44.748673 36.687615  0.95167893 1.4812747 44.748673 44.748673
#> 2 prior  2.010989  2.119266 -0.07443683 0.7480863  2.010989  2.010989
#>            smdw       vrw
#> 1 -8.388583e-16 1.9940894
#> 2 -3.052964e-16 0.6845163
```
