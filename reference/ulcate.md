# The U-learner of conditional average treatment effects

Estimates how the effect of a 0/1 treatment varies with the covariates
(Kunzel et al. 2019, sec. 2 and appendix): the outcome regression
`m(x) = E[y | x]` and the propensity `e(x) = P(d = 1 | x)` are learned
out of fold, as in
[`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
each fold predicted by learners fitted on the others; the pseudo-outcome
`u = (y - m) / (d - e)` has expectation `tau(x)` given `x` when the
treatment is unconfounded given `x`; and `tau` is the learner `lrnt`
fitted to `u` on all the rows. The mean of `tau` over the rows is
reported as `ate`.

## Usage

``` r
ulcate(
  y,
  d,
  x,
  fold,
  lrn = "lasso",
  lrnd = lrn,
  lrnt = "glm",
  clip = 0.01,
  alpha = 0.5,
  nlam = 100,
  nfin = 5,
  tol = 1e-07,
  scr = TRUE,
  hon = FALSE,
  ntree = 500,
  mtry = NULL,
  pred = NULL,
  stk = c("glm", "lasso", "forest"),
  nstk = 5,
  disc = FALSE,
  seed = 1
)

# S3 method for class 'ulcate'
print(x, ...)
```

## Arguments

- y:

  The outcome: one finite number per row of `x`.

- d:

  The treatment: 0 or 1, both present.

- x:

  The covariates: a numeric matrix, one row per observation.

- fold:

  The fold of each row, whole numbers from 1 to the number of folds (at
  least 2), each used.

- lrn:

  The learner of the outcome regression `m`: `"glm"`, `"lasso"`,
  `"ridge"`, `"enet"`, `"forest"` or `"stack"`, as in
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md).

- lrnd:

  The learner of the propensity `e`, from the same choices (a
  probability: logistic for `"glm"`).

- lrnt:

  The learner that fits `tau` to the pseudo-outcome, from the same
  choices; `"glm"` gives an effect linear in the covariates.

- clip:

  The smallest size of `d - e`: one number above 0 and below 0.5.

- alpha:

  The mixing parameter of `"enet"`, from 0 to 1.

- nlam:

  The number of penalties on each elastic net path, at most 67108864
  (2^26).

- nfin:

  The number of inner folds that choose the penalty, at least 2 and at
  most the rows of `x`.

- tol, scr:

  The convergence threshold and the screening of the elastic net
  learners, as in
  [`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md).

- hon:

  `TRUE` for honest forests, on half-samples without replacement with
  `imb = 0.05`, as 'grf' grows them (see
  [`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md)).

- ntree:

  The number of trees of `"forest"`: a whole number, at least 1.

- mtry:

  The columns tried at each split of `"forest"`: `NULL` for the default
  of
  [`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md),
  or a whole number (more than the columns means all of them).

- pred:

  `NULL` to learn `m` and `e`, or their predictions: a matrix with one
  row per observation and two columns, `m` then `e`.

- stk:

  The learners of `"stack"`, each named once, from the other choices of
  `lrn`.

- nstk:

  The number of inner folds of `"stack"`, at least 2.

- disc:

  `TRUE` for the discrete super learner in `"stack"` (the one learner of
  least cross-validated risk), `FALSE` for non-negative weights summing
  to 1.

- seed:

  The seed of the inner folds and of the forests: one number, at least
  0.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `ulcate`, a list: `tau`, the effect at each row;
`ate`, its mean; `u`, the pseudo-outcome; `pred`, the out-of-fold
predictions (columns `m` and `e`).

## Details

The denominator `d - e` is near 0 for a treated row with a propensity
near 1 or an untreated row with one near 0, where `u` explodes: it is
moved to `clip` (or `-clip`) when nearer 0 than that. The learners, the
folds, the clipping and the final fit are all compiled.

## References

Kunzel, S. R., Sekhon, J. S., Bickel, P. J. and Yu, B. (2019).
Metalearners for estimating heterogeneous treatment effects using
machine learning. *Proceedings of the National Academy of Sciences*,
116(10), 4156-4165.
[doi:10.1073/pnas.1804597116](https://doi.org/10.1073/pnas.1804597116)

Nie, X. and Wager, S. (2021). Quasi-oracle estimation of heterogeneous
treatment effects. *Biometrika*, 108(2), 299-319.
[doi:10.1093/biomet/asaa076](https://doi.org/10.1093/biomet/asaa076)

## Examples

``` r
set.seed(6)
n <- 400
x <- matrix(rnorm(n * 2), n)
d <- rbinom(n, 1, plogis(0.5 * x[, 1]))
y <- x[, 2] + d * (1 + x[, 1]) + rnorm(n)
fold <- sample(rep(1:4, length.out = n))
r <- ulcate(y, d, x, fold, lrn = "glm", lrnd = "glm")
r
#> ulcate: U-learner of the effect of d, 400 rows; mean effect 0.9368
#> tau from -2.416 to 3.844
```
