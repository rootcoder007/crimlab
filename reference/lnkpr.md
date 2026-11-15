# Crime linkage: which crimes belong to one series

Asks whether two crimes were committed by the same offender, from how
far apart they are in space and time and how alike the offender behaved
(the modus operandi: entry method, property taken, a weapon shown). For
every pair of crimes it gives the distance, the time between them and
the Jaccard similarity of their behaviours (the behaviours present in
both over those present in either; 0 when neither shows any). On the
pairs of solved crimes, linked when they belong to the same series, a
logistic regression of being linked on the three measures is fitted and
each measure's power to tell linked from unlinked pairs is the area
under its ROC curve (AUC: the chance that a linked pair is closer,
sooner or more alike than an unlinked one; 0.5 is no better than
chance), as Bennell and Canter (2002) did for commercial burglaries. The
fitted model then scores every pair, including those with an unsolved
crime: `prob` is the estimated chance that the two are linked.

## Usage

``` r
lnkpr(x, y = NULL, t, series, mo = NULL)
```

## Arguments

- x, y:

  The points in planar coordinates (project longitude and latitude with
  [`eqproj()`](https://rootcoder007.github.io/crimlab/reference/eqproj.md)
  first): `x` and `y` numeric vectors of equal length, or `x` a matrix
  or data frame with columns named `x` and `y`. At least two points; no
  missing values.

- t:

  The time of each crime (days, say): finite numbers.

- series:

  The series of each solved crime: labels, `NA` for an unsolved crime.

- mo:

  `NULL` or the behaviours: a 0/1 or logical matrix, or a data frame or
  list of 0/1 columns, one row per crime.

## Value

An object of class `lnkpr`, a list: `coef`, the logistic regression
(`term`, `est`, `se`, `z`, `p`); `auc`, a named vector (`dist`, `dt`,
`jac`, `model`); `pairs`, one row per pair of crimes (`i`, `j`, their
rows, `dist`, `dt`, `jac`, `linked`, `NA` when either is unsolved, and
`prob`); `dev`, the model's deviance; and `nmo`, the number of
behaviours.

## Details

The probabilities are those of the pairs in the data: series with many
crimes give many linked pairs, and the share linked among pairs depends
on how many series there are. The AUC of the model is in-sample and so
optimistic.

## References

Bennell, C. and Canter, D. V. (2002). Linking commercial burglaries by
modus operandi: tests using regression and ROC analysis. *Science &
Justice*, 42(3), 153-164.
[doi:10.1016/S1355-0306(02)71820-0](https://doi.org/10.1016/S1355-0306%2802%2971820-0)

## Examples

``` r
set.seed(9)
s <- rep(1:6, each = 5)
cx <- runif(6, 0, 5000)[s] + rnorm(30, 0, 1500)
cy <- runif(6, 0, 5000)[s] + rnorm(30, 0, 1500)
day <- runif(6, 0, 300)[s] + rnorm(30, 0, 60)
mo <- matrix(rbinom(30 * 5, 1, 0.3), 30)
sr <- ifelse(runif(30) < 0.2, NA, s)
r <- lnkpr(cbind(x = cx, y = cy), t = day, series = sr, mo = mo)
r$coef
#>          term           est           se          z            p
#> 1 (Intercept)  0.1466816977 0.4786319854  0.3064603 7.592542e-01
#> 2        dist -0.0001929138 0.0001228189 -1.5707178 1.162482e-01
#> 3          dt -0.0157565696 0.0033891231 -4.6491582 3.332925e-06
#> 4         jac  0.5344187921 0.8865379558  0.6028155 5.466315e-01
r$auc
#>      dist        dt       jac     model 
#> 0.5950469 0.7775642 0.5022617 0.7888726 
```
