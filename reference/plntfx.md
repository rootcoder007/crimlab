# Does the analysis recover a planted effect?

Plants effects of known size in the real data and checks that the
analysis finds them: the outcome is shuffled across rows, which removes
whatever association it had with the treatment and the covariates while
keeping its distribution, an effect `eff` is added to the treated
(`y* = shuffled y + eff d`), and the analysis is run again, `nrep` times
for each effect. An analysis that is biased, or whose intervals are too
narrow, shows it here on data where the truth is known, a positive
control in the sense of Lipsitch, Tchetgen Tchetgen and Cohen (2010).

## Usage

``` r
plntfx(
  data,
  y,
  d,
  x = NULL,
  eff,
  nrep = 200,
  se = c("hc1", "iid"),
  level = 0.95,
  tol = NULL,
  seed = NULL
)

# S3 method for class 'plntfx'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit.

- y:

  The name of the outcome column (numbers).

- d:

  The name of the treatment column (numbers; it must vary).

- x:

  `NULL` or the names of the covariate columns (numbers); for
  [`print()`](https://rdrr.io/r/base/print.html), an object of class
  `plntfx`.

- eff:

  The effects to plant, in the outcome's units: finite numbers.

- nrep:

  The shuffles per effect: a whole number from 2 to 1000000.

- se:

  `"hc1"` or `"iid"`.

- level:

  The level of the intervals: one number strictly between 0 and 1.

- tol:

  `NULL` (three Monte Carlo standard errors) or the largest `|bias|`
  accepted: one number from 0.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `plntfx`, a list: `tab`, a data frame with one row
per effect (`eff`, `mean`, `bias`, `sd`, `rmse`, `cover`, `rej`, `tol`,
`ok`); and `nrep`, `se`, `level` and `seed`.

## Details

The analysis is the coefficient of `d` in the least squares fit of the
outcome on an intercept, `x` and `d`, with a classical or HC1 standard
error and Student's t interval at `level`, as
[`falsctl()`](https://rootcoder007.github.io/crimlab/reference/falsctl.md)
fits it. For each effect: `mean` and `sd` of the estimates,
`bias = mean - eff`, `rmse` (about `eff`), `cover` (the share of
intervals holding `eff`), `rej` (the share rejecting 0 at `1 - level`),
and `ok`, whether `|bias|` is at most `tol`. By default `tol` is three
Monte Carlo standard errors of the mean, `3 sd / sqrt(nrep)`. One
shuffle serves every effect, and the shuffles come from the package's
generator.

## References

Lipsitch, M., Tchetgen Tchetgen, E. and Cohen, T. (2010). Negative
controls: a tool for detecting confounding and bias in observational
studies. *Epidemiology*, 21(3), 383-388.
[doi:10.1097/EDE.0b013e3181d61eeb](https://doi.org/10.1097/EDE.0b013e3181d61eeb)

Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine*,
38(11), 2074-2102.
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086)

## Examples

``` r
set.seed(6)
d <- data.frame(age = runif(200, 18, 60), t = rbinom(200, 1, 0.4))
d$y <- 0.02 * d$age + rnorm(200)
plntfx(d, "y", "t", x = "age", eff = c(0, 0.25, 0.5), nrep = 200,
       seed = 1)
#> plntfx: 3 of 3 planted effect(s) recovered within tolerance over 200 shuffles
#>   eff      mean      bias     sd   rmse cover   rej     tol   ok
#>  0.00 -0.002595 -0.002595 0.1514 0.1511  0.93 0.070 0.03212 TRUE
#>  0.25  0.247405 -0.002595 0.1514 0.1511  0.93 0.395 0.03212 TRUE
#>  0.50  0.497405 -0.002595 0.1514 0.1511  0.93 0.920 0.03212 TRUE
```
