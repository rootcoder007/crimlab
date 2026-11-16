# Falsification with negative-control outcomes

Estimates the "effect" of the treatment on outcomes it cannot have
affected (negative controls: an outcome measured before treatment, or
one with no plausible path from it), by the same regression as the
outcome of interest. An association with a negative control points to
confounding or another bias that the adjustment did not remove
(Lipsitch, Tchetgen Tchetgen and Cohen 2010; Shi, Miao and Tchetgen
Tchetgen 2020).

## Usage

``` r
falsctl(
  data,
  y = NULL,
  d,
  nc,
  x = NULL,
  se = c("hc1", "iid"),
  nperm = 0,
  adj = c("holm", "bonferroni", "hochberg", "BH", "BY", "none"),
  alpha = 0.05,
  level = 0.95,
  seed = NULL
)

# S3 method for class 'falsctl'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  `NULL` or the name of the outcome of interest (numbers).

- d:

  The name of the treatment column (numbers; it must vary).

- nc:

  The names of the negative-control outcome columns (numbers).

- x:

  `NULL` or the names of the covariate columns (numbers); for
  [`print()`](https://rdrr.io/r/base/print.html), an object of class
  `falsctl`.

- se:

  `"hc1"` or `"iid"`.

- nperm:

  The permutations: a whole number from 0 (none) to 1000000.

- adj:

  The adjustment of the controls' p-values: `"holm"`, `"bonferroni"`,
  `"hochberg"`, `"BH"`, `"BY"` or `"none"`.

- alpha:

  The level at which a control is flagged: one number strictly between 0
  and 1.

- level:

  The level of the intervals: one number strictly between 0 and 1.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647; not used without permutations.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `falsctl`, a list: `tab`, a data frame with one row
per outcome (`outcome`, `role`, `primary` or `control`, `est`, `se`,
`stat`, `p`, `lo`, `hi`, `pperm`, `padj`; `NA` where not computed);
`nflag`; `ncontrol`; and `se`, `adj`, `alpha`, `level`, `nperm`, `seed`
(`NA` without permutations).

## Details

Each outcome is regressed on an intercept, the covariates `x` and the
treatment `d` by least squares; `est` is the coefficient of `d` with a
classical (`se = "iid"`) or HC1 (`"hc1"`) standard error, Student's t on
the residual degrees of freedom and its interval at `level`. With
`nperm > 0` each outcome also has the permutation p-value of Kennedy
(1995): the residuals of the outcome on the covariates are permuted
against those of `d`, one permutation of the rows for every outcome,
`pperm = (1 + k) / (1 + nperm)`. The controls' p-values (the permutation
ones when drawn) are adjusted for their number by `adj` (the methods of
[`p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)); `nflag` counts
the controls whose adjusted p-value is at most `alpha`.

## References

Lipsitch, M., Tchetgen Tchetgen, E. and Cohen, T. (2010). Negative
controls: a tool for detecting confounding and bias in observational
studies. *Epidemiology*, 21(3), 383-388.
[doi:10.1097/EDE.0b013e3181d61eeb](https://doi.org/10.1097/EDE.0b013e3181d61eeb)

Shi, X., Miao, W. and Tchetgen Tchetgen, E. (2020). A selective review
of negative control methods in epidemiology. *Current Epidemiology
Reports*, 7(4), 190-202.
[doi:10.1007/s40471-020-00243-4](https://doi.org/10.1007/s40471-020-00243-4)

Kennedy, P. E. (1995). Randomization tests in econometrics. *Journal of
Business and Economic Statistics*, 13(1), 85-94.
[doi:10.1080/07350015.1995.10524581](https://doi.org/10.1080/07350015.1995.10524581)

## Examples

``` r
set.seed(5)
n <- 300
d <- data.frame(age = runif(n, 18, 60), u = rnorm(n))
d$t <- rbinom(n, 1, plogis(0.8 * d$u))
d$arrest <- 0.4 * d$t + 0.5 * d$u + rnorm(n)
d$prior <- 0.5 * d$u + rnorm(n)     # before treatment: a control
d$height <- rnorm(n)                # no path from treatment
falsctl(d, "arrest", "t", c("prior", "height"), x = "age")
#> falsctl: 1 of 2 negative control(s) associated with the treatment at alpha 0.05 (holm)
#>  outcome    role     est     se   stat         p      lo       hi pperm
#>   arrest primary  0.8024 0.1308  6.134 2.730e-09  0.5450 1.059832    NA
#>    prior control  0.3816 0.1282  2.976 3.159e-03  0.1293 0.633989    NA
#>   height control -0.2236 0.1171 -1.909 5.720e-02 -0.4541 0.006887    NA
#>      padj
#>        NA
#>  0.006318
#>  0.057202
```
