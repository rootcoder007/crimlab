# Random forest for regression and probabilities

Grows a random forest (Breiman 2001) on `x` and `y` and returns its mean
prediction at the rows of `newx` and its out-of-bag prediction at each
row of `x`. Each tree is grown on a bag of rows drawn from `x`, as deep
as `minnode` allows, choosing at each node the best of `mtry` columns
drawn at random by the reduction in the sum of squares; for a 0/1
outcome this is the Gini split, and the forest's mean is a probability,
as in a probability forest of 'ranger'.

## Usage

``` r
rfor(
  x,
  y,
  newx = x,
  ntree = 500,
  mtry = NULL,
  minnode = 5,
  frac = 1,
  repl = TRUE,
  hon = FALSE,
  imb = 0,
  seed = 1
)
```

## Arguments

- x:

  The predictors: a numeric matrix with at least two rows, no missing
  values.

- y:

  The outcome: one finite number per row of `x`; 0 or 1 for a
  probability.

- newx:

  The rows to predict: a numeric matrix with the columns of `x`.

- ntree:

  The number of trees: one whole number, at least 1.

- mtry:

  The columns tried at each split: `NULL` for
  `min(ceiling( sqrt(p) + 20), p)`, `p` the number of columns (the rule
  of 'grf'), or a whole number from 1 to the number of columns.

- minnode:

  Nodes of at most this many rows are not split: one whole number, at
  least 1 (5 is the default of 'ranger' for regression, 10 for
  probabilities).

- frac:

  The size of each bag as a share of the rows of `x`: above 0, at most
  1.

- repl:

  `TRUE` to draw bags with replacement, `FALSE` without.

- hon:

  `TRUE` for honest trees (Wager and Athey 2018): each bag is shuffled,
  its first half places the splits and its second half sets the leaf
  means (a leaf without such rows takes the mean of its deepest ancestor
  that has some). Needs bags of at least two rows.

- imb:

  Each child of a split must hold at least this share of its node's rows
  (and one row): one number from 0 to below 0.5 ('grf' uses 0.05).

- seed:

  The seed of the draws: one number, at least 0. The same seed and data
  give the same forest.

## Value

A list: `pred`, the predictions at the rows of `newx`, and `oob`, the
out-of-bag prediction at each row of `x` (`NaN` for a row that every
tree's bag holds).

## Details

Each column is cut once, at its own order statistics, into at most 256
bins: every split between neighbouring values when a column has at most
256 distinct values (as an exact search finds), and splits at 255
quantiles otherwise. A node's best split then comes from one pass over
its rows, and the forest depends on each column only through the order
of its values: a strictly increasing transform of a column changes
nothing. The forest is not bit-for-bit that of 'ranger' (it draws its
own random numbers and uses quantile bins on columns of many values);
tests compare their errors on known functions instead.

## Defaults

The defaults were chosen by test error against the true function over
four data-generating functions (a sparse linear one, Friedman's, a
step-and-curve one and a probability), p = 20 and 200, ten seeds each
(`dev/rftune.R`). Like for like (the same trees, `mtry`, node size,
half-samples without replacement, `imb = 0.05`, honest or not), `rfor()`
and 'grf' have test errors within 1 to 3 per cent of each other, inside
two standard errors in almost every cell. The one setting that mattered
was `mtry`: the rule of 'grf' lowers the error against the square root
of `p` by 30 to 50 per cent on the sparse and step functions, so it is
the default. Bags of all rows drawn with replacement gave the lowest
error with that `mtry` in most cells; honest trees and half-samples
raised the error of the predictions (honesty serves inference, not
prediction), so both are options, off by default. Nodes of 5 rows are
left unsplit, as both packages do.

## References

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
[doi:10.1023/A:1010933404324](https://doi.org/10.1023/A%3A1010933404324)

Wager, S. and Athey, S. (2018). Estimation and inference of
heterogeneous treatment effects using random forests. *Journal of the
American Statistical Association*, 113(523), 1228-1242.
[doi:10.1080/01621459.2017.1319839](https://doi.org/10.1080/01621459.2017.1319839)

## Examples

``` r
set.seed(1)
x <- matrix(rnorm(400 * 4), 400)
y <- sin(2 * x[, 1]) + x[, 2] + rnorm(400, sd = 0.3)
f <- rfor(x, y, ntree = 100)
mean((f$oob - y)^2)
#> [1] 0.184494
```
