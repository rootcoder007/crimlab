# Generalisability weights for a trial

The weights that make the units of a randomised trial stand for a target
population, and how far the trial is from it. The data hold the trial's
units (`s = 1`) and a sample of the target (`s = 0`) with the same
covariates. The sampling score `ps`, the chance of being in the trial
given the covariates, is a logistic regression of `s` on `x`, trimmed to
`[trim, 1 - trim]`.

## Usage

``` r
gensw(data, s, x, target = c("generalize", "transport"), trim = 0, bins = 10)

# S3 method for class 'gensw'
print(x, ...)
```

## Arguments

- data:

  A data frame, one row per unit of the trial or the target.

- s:

  The name of the column marking the trial's units: 1 in the trial, 0
  outside it (or `TRUE`, `FALSE`), both present.

- x:

  The names of the covariate columns (numbers) of the sampling model;
  for [`print()`](https://rdrr.io/r/base/print.html), an object of class
  `gensw`.

- target:

  `"generalize"` or `"transport"`.

- trim:

  The bound of the scores, from 0 to below 0.5.

- bins:

  The bins of Tipton's index: a whole number from 1 to 10000.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `gensw`, a list: `ps` and `w`, one per row; `coef`,
the sampling model's coefficients; `diag`, a named vector (`ntrial`,
`ntarget`, `smd`, `dmean`, `gindex`, `ess`); `bal`, a data frame (`var`,
`target`, `trial`, `wtrial`, `smd`, `wsmd`); and `target`.

## Details

With `target = "generalize"` the target is the population the trial was
drawn from, all the rows, and a trial unit's weight is `1 / ps`; with
`"transport"` it is the units outside the trial, and the weight is the
odds `(1 - ps) / ps` (Stuart et al. 2011; Westreich et al. 2017; Tipton
and Hartman 2023, eq. 3.11-3.12). Rows outside the trial have weight 0.

The diagnostics: `smd`, the difference of the mean logits of the scores
between the trial and the target over their standard deviation in the
target (Stuart et al. 2011); `dmean`, the difference of the mean scores;
`gindex`, Tipton's (2014) generalisability index, the sum over `bins`
equal-width bins of the scores of the square root of the trial's share
times the target's, 1 when they are spread alike and 0 when they never
meet (above 0.9 is very high, below 0.5 low); `ess`, Kish's effective
number of trial units under the weights. `bal` gives each covariate's
mean in the target and in the trial, unweighted and weighted, with the
standardised differences from the target (over its standard deviation in
the target).

## References

Stuart, E. A., Cole, S. R., Bradshaw, C. P. and Leaf, P. J. (2011). The
use of propensity scores to assess the generalizability of results from
randomized trials. *Journal of the Royal Statistical Society, Series A*,
174(2), 369-386.
[doi:10.1111/j.1467-985X.2010.00673.x](https://doi.org/10.1111/j.1467-985X.2010.00673.x)

Tipton, E. (2014). How generalizable is your experiment? An index for
comparing experimental samples and populations. *Journal of Educational
and Behavioral Statistics*, 39(6), 478-501.
[doi:10.3102/1076998614558486](https://doi.org/10.3102/1076998614558486)

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
set.seed(8)
n <- 600
d <- data.frame(age = rnorm(n, 35, 10), prior = rpois(n, 1.5))
d$s <- rbinom(n, 1, plogis(-1 + 0.04 * (d$age - 35) + 0.3 * d$prior))
g <- gensw(d, "s", c("age", "prior"), target = "transport")
g
#> gensw: 223 trial units to transport to 377; logit SMD 0.47, Tipton's index 0.96, effective size 183.7
#>    var target  trial wtrial    smd      wsmd
#>    age 33.012 36.075 32.979 0.3000 -0.003267
#>  prior  1.358  1.771  1.317 0.3693 -0.036677
g$bal
#>     var   target    trial    wtrial       smd         wsmd
#> 1   age 33.01239 36.07538 32.979032 0.2999789 -0.003266946
#> 2 prior  1.35809  1.77130  1.317055 0.3693281 -0.036677349
```
