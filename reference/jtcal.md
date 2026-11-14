# Empirical distance decay from solved trips

A kernel density estimate of trip lengths from known offender trips, to
calibrate a decay function (Levine 2009; Silverman 1986). With
`refl = TRUE` the kernels are reflected at 0, so no mass falls below 0
and the estimate integrates to 1 over `[0, Inf)`. Tabulated at evenly
spaced distances it is a `"table"` decay for
[`jtdecay()`](https://rootcoder007.github.io/crimlab/reference/jtdecay.md)
and
[`jtsurf()`](https://rootcoder007.github.io/crimlab/reference/jtsurf.md).

## Usage

``` r
jtcal(d, t, h, w = NULL, kern = "normal", refl = TRUE)
```

## Arguments

- d:

  Trip lengths: numbers from 0.

- t:

  The distances to evaluate at (from 0 when `refl = TRUE`).

- h:

  The bandwidth: one positive number.

- w:

  `NULL` for equal weights, or one weight from 0 per trip, with a
  positive sum.

- kern:

  The kernel: `"normal"`, `"quartic"`, `"triangular"`, `"uniform"` or
  `"laplace"` (double exponential).

- refl:

  `TRUE` to reflect at 0.

## Value

A numeric vector of densities, one per value of `t`.

## References

Silverman, B. W. (1986). *Density Estimation for Statistics and Data
Analysis*. Chapman and Hall; reissued 2018 by Routledge.
[doi:10.1201/9781315140919](https://doi.org/10.1201/9781315140919)

Levine, N. (2009). Introduction to the special issue on Bayesian
journey-to-crime modelling. *Journal of Investigative Psychology and
Offender Profiling*, 6(3), 167-185.
[doi:10.1002/jip.107](https://doi.org/10.1002/jip.107)

## Examples

``` r
set.seed(1)
trips <- rexp(200, 1 / 2)
jtcal(trips, t = seq(0, 10, 0.5), h = 0.5)
#>  [1] 0.341941995 0.334326174 0.303888364 0.263429537 0.223093541 0.178452942
#>  [7] 0.129191709 0.095269714 0.079083025 0.063478487 0.046114240 0.036675173
#> [13] 0.025708978 0.013938566 0.006084284 0.003949455 0.004978870 0.005423329
#> [19] 0.005837164 0.005507356 0.003471230
```
