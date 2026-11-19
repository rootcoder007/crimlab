# The parts of a fitted regression

The fits of
[`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md),
[`expgmm()`](https://rootcoder007.github.io/crimlab/reference/expgmm.md),
[`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md),
[`gmmpnl()`](https://rootcoder007.github.io/crimlab/reference/gmmpnl.md),
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md),
[`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)
and
[`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
are lists, and the usual accessors read them:
[`coef()`](https://rdrr.io/r/stats/coef.html) gives the estimates by
term, [`confint()`](https://rdrr.io/r/stats/confint.html) their
intervals, [`vcov()`](https://rdrr.io/r/stats/vcov.html) the covariance
matrix of the estimates and
[`nobs()`](https://rdrr.io/r/stats/nobs.html) the rows the fit used.

## Usage

``` r
# S3 method for class 'ppmlfe'
coef(object, ...)

# S3 method for class 'ivfe'
coef(object, ...)

# S3 method for class 'ppmlcf'
coef(object, ...)

# S3 method for class 'expgmm'
coef(object, ...)

# S3 method for class 'fepnl'
coef(object, ...)

# S3 method for class 'gmmpnl'
coef(object, ...)

# S3 method for class 'didfe'
coef(object, ...)

# S3 method for class 'rdloc'
coef(object, ...)

# S3 method for class 'sarml'
coef(object, ...)

# S3 method for class 'ppmlfe'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'ivfe'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'ppmlcf'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'expgmm'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'fepnl'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'gmmpnl'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'didfe'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'rdloc'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'sarml'
confint(object, parm = NULL, level = NULL, ...)

# S3 method for class 'ppmlfe'
vcov(object, ...)

# S3 method for class 'ivfe'
vcov(object, ...)

# S3 method for class 'expgmm'
vcov(object, ...)

# S3 method for class 'fepnl'
vcov(object, ...)

# S3 method for class 'gmmpnl'
vcov(object, ...)

# S3 method for class 'didfe'
vcov(object, ...)

# S3 method for class 'sarml'
vcov(object, ...)

# S3 method for class 'ppmlfe'
nobs(object, ...)

# S3 method for class 'ivfe'
nobs(object, ...)

# S3 method for class 'ppmlcf'
nobs(object, ...)

# S3 method for class 'expgmm'
nobs(object, ...)

# S3 method for class 'fepnl'
nobs(object, ...)

# S3 method for class 'gmmpnl'
nobs(object, ...)

# S3 method for class 'didfe'
nobs(object, ...)

# S3 method for class 'rdloc'
nobs(object, ...)

# S3 method for class 'sarml'
nobs(object, ...)
```

## Arguments

- object:

  A fit of one of the functions above.

- ...:

  Nothing; the generics of stats have it.

- parm:

  `NULL` (every term) or the terms wanted: names of terms or their
  positions.

- level:

  The level of the intervals: `NULL` (the fit's own; 0.95 for
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md))
  or one number between 0 and 1. For every fit but that of
  [`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
  it must be the fit's level.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): a named numeric vector,
one estimate per term (per method for
[`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)).
[`confint()`](https://rdrr.io/r/stats/confint.html): a matrix with a row
per term and the columns `lower` and `upper`, with the level as the
attribute `"level"`. [`vcov()`](https://rdrr.io/r/stats/vcov.html): a
square matrix named by the terms.
[`nobs()`](https://rdrr.io/r/stats/nobs.html): one whole number.

## Details

The intervals are those the fit computed, at the fit's own level (the
bias-corrected interval of
[`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md),
the bootstrap interval of
[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)):
they are read from the fit, not rebuilt from the standard errors, so
[`confint()`](https://rdrr.io/r/stats/confint.html) takes no other
level; fit again with the `level` wanted.
[`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
stores no interval, and its
[`confint()`](https://rdrr.io/r/stats/confint.html) is the normal one,
estimate plus or minus the quantile times the standard error, at any
`level`.

[`ppmlcf()`](https://rootcoder007.github.io/crimlab/reference/ppmlcf.md)
(bootstrap) and
[`rdloc()`](https://rootcoder007.github.io/crimlab/reference/rdloc.md)
(one estimate by two methods) store no covariance matrix, so they have
no [`vcov()`](https://rdrr.io/r/stats/vcov.html). The matrix of
[`sarml()`](https://rootcoder007.github.io/crimlab/reference/sarml.md)
has a row for `rho` and one for `lambda` whichever the model holds
(`NaN` for the one it does not), as its help page says.

## How a model is given

By the names of columns of `data`, or by vectors, each in the argument
of its role: the outcome, the regressors, the instruments, the fixed
effects, the clusters. There is no formula interface: a formula has no
place for those roles, and the fixed effects are absorbed, not expanded
into columns of a model matrix. Regressors are numeric columns; a column
of labels is taken only as a fixed effect or a cluster. Missing and
infinite values are refused, each in its own words, and no option drops
them: the caller subsets first.

## What a fit keeps

The terms keep the names of the input columns, in
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html). Row names are not kept: a
fit has one row per term, not per observation. Fitted values need the
fixed effects, which are absorbed and never estimated, so the fits do
not return them;
[`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md)
and
[`gmmpnl()`](https://rootcoder007.github.io/crimlab/reference/gmmpnl.md)
return their residuals (`resid`). Collinear regressors, and regressors
the fixed effects absorb, stop the fit with a message that says so. A
count fit that has not converged warns, and records it (`conv`, `iter`);
`tol` and `maxit` set the stopping rule.

## Speed

The time of a fit grows in proportion to the number of rows: ten times
the rows take about ten times as long (measured for
[`ivfe()`](https://rootcoder007.github.io/crimlab/reference/ivfe.md),
[`ppmlfe()`](https://rootcoder007.github.io/crimlab/reference/ppmlfe.md),
[`fepnl()`](https://rootcoder007.github.io/crimlab/reference/fepnl.md)
and
[`didfe()`](https://rootcoder007.github.io/crimlab/reference/didfe.md)
from a thousand to a hundred thousand rows).

## Examples

``` r
set.seed(1)
d <- data.frame(u = rep(1:20, each = 6), x = rnorm(120))
d$y <- 1 + 0.5 * d$x + rnorm(120)
f <- ivfe(d, "y", x = "x", fe = "u")
coef(f)
#>         x 
#> 0.4803809 
confint(f)
#>       lower     upper
#> x 0.1880123 0.7727495
#> attr(,"level")
#> [1] 0.95
vcov(f)
#>            x
#> x 0.01951251
nobs(f)
#> [1] 120
```
