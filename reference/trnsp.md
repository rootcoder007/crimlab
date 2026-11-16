# Transport a trial's effect to a target population

The average treatment effect of a randomised trial carried over to a
target population whose covariates differ from the trial's, by the
outcome model and by inverse odds weighting (Dahabreh et al. 2020;
Westreich et al. 2017). The data hold the trial's units (`s = 1`), with
outcome and treatment, and a sample of the target (`s = 0`), with the
covariates only; `y` and `d` are not read outside the trial and may be
`NA` there.

## Usage

``` r
trnsp(
  data,
  y,
  d,
  s,
  x,
  target = c("transport", "generalize"),
  trim = 0,
  nboot = 200,
  level = 0.95,
  seed = NULL
)

# S3 method for class 'trnsp'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit of the trial or the target.

- y:

  The name of the outcome column (numbers; read in the trial).

- d:

  The name of the treatment column: 0 or 1 in the trial, both present.

- s:

  The name of the column marking the trial's units: 1 in the trial, 0
  outside it (or `TRUE`, `FALSE`), both present.

- x:

  The names of the covariate columns (numbers), complete in every row;
  for [`print()`](https://rdrr.io/r/base/print.html), an object of class
  `trnsp`.

- target:

  `"transport"` or `"generalize"`.

- trim:

  The bound of the sampling scores, from 0 to below 0.5.

- nboot:

  The bootstrap draws: a whole number from 0 (no standard errors) to
  100000.

- level:

  The level of the intervals: one number strictly between 0 and 1.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647; not used without draws.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `trnsp`, a list: `tab`, a data frame with rows `om`
and `iow` (`method`, `est`, `se`, `lo`, `hi`); `ps` and `w`, the
sampling scores and weights of every row; `draws` (a matrix, `NA` for a
draw left out); `nfail`; and `target`, `nboot`, `level`, `seed`.

## Details

The target is the rows outside the trial (`target = "transport"`) or all
rows, trial included (`"generalize"`). The outcome model fits the
outcome on an intercept and `x` by least squares in each arm of the
trial and averages the difference of the two predictions over the target
(`om`). Inverse odds weighting gives each trial unit the weight of
[`gensw()`](https://rootcoder007.github.io/crimlab/reference/gensw.md)
(`(1 - ps) / ps` to transport, `1 / ps` to generalise, the sampling
score `ps` trimmed to `[trim, 1 - trim]`) and takes the difference of
the arms' weighted means, the weights normalised within each arm (`iow`;
Tipton and Hartman 2023, eq. 3.10). The first is right when the outcome
model is, the second when the sampling model is. Standard errors come
from `nboot` bootstrap draws of the rows within the trial and within the
rest, a draw whose fit fails left out (`nfail`), with normal intervals
at `level`.

## References

Dahabreh, I. J., Robertson, S. E., Steingrimsson, J. A., Stuart, E. A.
and Hernan, M. A. (2020). Extending inferences from a randomized trial
to a new target population. *Statistics in Medicine*, 39(14), 1999-2014.
[doi:10.1002/sim.8426](https://doi.org/10.1002/sim.8426)

Westreich, D., Edwards, J. K., Lesko, C. R., Stuart, E. and Cole, S. R.
(2017). Transportability of trial results using inverse odds of sampling
weights. *American Journal of Epidemiology*, 186(8), 1010-1014.
[doi:10.1093/aje/kwx164](https://doi.org/10.1093/aje/kwx164)

Tipton, E. and Hartman, E. (2023). Generalizability and
transportability. In Zubizarreta, J. R., Stuart, E. A., Small, D. S. and
Rosenbaum, P. R. (eds), *Handbook of Matching and Weighting Adjustments
for Causal Inference*, ch. 3. Chapman and Hall/CRC.
[doi:10.1201/9781003102670](https://doi.org/10.1201/9781003102670)

## Examples

``` r
set.seed(9)
n <- 600
d <- data.frame(age = rnorm(n, 35, 10))
d$s <- rbinom(n, 1, plogis(-0.05 * (d$age - 35)))
d$t <- ifelse(d$s == 1, rbinom(n, 1, 0.5), NA)
d$y <- ifelse(d$s == 1, 0.5 * d$t + 0.02 * d$t * d$age + rnorm(n), NA)
trnsp(d, "y", "t", "s", "age", nboot = 100, seed = 1)
#> trnsp: the trial's effect transported to the rows outside it, 100 bootstrap draws (0 failed)
#>  method   est     se     lo    hi
#>      om 1.182 0.1282 0.9308 1.433
#>     iow 1.197 0.1288 0.9442 1.449
```
