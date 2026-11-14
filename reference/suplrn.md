# Stacking: the super learner

Combines the learners of the package (least squares or logistic
regression, the elastic net of
[`enpath()`](https://rootcoder007.github.io/crimlab/reference/enpath.md),
the random forest of
[`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md))
into one predictor by stacking (van der Laan, Polley and Hubbard 2007).
Each learner is fitted K times, on the rows outside each fold, and
predicts the fold left out; the meta-weights are the non-negative least
squares fit of `y` on these cross-fitted predictions, scaled to sum to 1
(the `method.NNLS` of 'SuperLearner'). With `disc = TRUE` the stack is
the discrete super learner: the one learner of least cross-validated
risk. The predictions are then the weighted sum of the learners fitted
on all the rows. Everything, the folds included, is computed in the
compiled core and is the same for the same `seed`; with one learner the
stack is that learner, bit for bit.

## Usage

``` r
suplrn(
  x,
  y,
  fam = "gaussian",
  lrns = c("glm", "lasso", "forest"),
  nfold = 5,
  fold = NULL,
  disc = FALSE,
  newx = NULL,
  alpha = 0.5,
  nlam = 100,
  nfin = 5,
  ntree = 500,
  mtry = NULL,
  seed = 1
)
```

## Arguments

- x:

  The predictors: a numeric matrix with at least four rows, no missing
  values.

- y:

  The outcome: one finite number per row of `x`; 0 or 1 for the binomial
  family.

- fam:

  `"gaussian"` or `"binomial"`.

- lrns:

  The learners to stack, each named once: from `"glm"`, `"lasso"`,
  `"ridge"`, `"enet"` and `"forest"`, as in
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md).

- nfold:

  The number of folds when `fold` is `NULL`: one whole number from 2 to
  the rows of `x`.

- fold:

  `NULL` to draw the folds from `seed` (a shuffle of the rows dealt
  round robin, so that fold sizes differ by at most one), or the fold of
  each row: whole numbers from 1 to the number of folds, each used.

- disc:

  `TRUE` for the discrete super learner (all the weight on the learner
  of least risk), `FALSE` for the non-negative weights.

- newx:

  `NULL` to predict the rows of `x`, or the rows to predict: a numeric
  matrix with the columns of `x`.

- alpha, nlam:

  The mixing parameter of `"enet"`, from 0 to 1, and the number of
  penalties on each elastic net path, as in
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md).

- nfin:

  The number of inner folds that choose the penalty of the elastic net,
  at least 2.

- ntree:

  The number of trees of `"forest"`: a whole number, at least 1.

- mtry:

  The columns tried at each split of `"forest"`: `NULL` for the default
  of
  [`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md),
  or a whole number.

- seed:

  The seed of the folds, the inner folds and the forests: one number, at
  least 0.

## Value

A list: `wt`, the meta-weights (non-negative, summing to 1) and `risk`,
the cross-validated mean squared error of each learner, both named by
learner; `fold`, the fold of each row; `oof`, the cross-fitted
predictions, one column per learner; and `pred`, the stack's predictions
at the rows of `newx`.

## Details

For a 0/1 outcome the learners predict probabilities (logistic
regression, the binomial elastic net, the forest's mean), the risk is
the Brier score and the weighted sum is again a probability.

[`dmlx()`](https://rootcoder007.github.io/crimlab/reference/dmlx.md)
takes the same stack as a nuisance learner (`lrn = "stack"`).

## References

van der Laan, M. J., Polley, E. C. and Hubbard, A. E. (2007). Super
learner. *Statistical Applications in Genetics and Molecular Biology*,
6(1), 25.
[doi:10.2202/1544-6115.1309](https://doi.org/10.2202/1544-6115.1309)

Lawson, C. L. and Hanson, R. J. (1974). *Solving Least Squares
Problems*. Prentice-Hall, Englewood Cliffs, NJ.

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(300 * 4), 300)
y <- sin(2 * x[, 1]) + x[, 2] + rnorm(300, sd = 0.3)
s <- suplrn(x, y, lrns = c("glm", "lasso", "forest"), ntree = 50)
s$wt
#>       glm     lasso    forest 
#> 0.0000000 0.0234579 0.9765421 
s$risk
#>       glm     lasso    forest 
#> 0.5327638 0.5345319 0.2611013 
```
