# Permutation importance of a prediction function

How much worse a risk tool predicts when one of its inputs is made
useless: the column is shuffled across the rows, which keeps its
distribution but breaks its link with the outcome and with the other
columns, and the loss of the predictions against `y` is measured again
(Breiman 2001; Fisher, Rudin and Dominici 2019, who call the ratio the
model reliance). Each column is shuffled `nrep` times and the mean and
the standard deviation of the increases are reported.

## Usage

``` r
xpimp(
  x,
  f,
  y,
  nrep = 10,
  loss = c("se", "ae"),
  ratio = FALSE,
  seed = 1,
  lctl = list()
)

# S3 method for class 'xpimp'
print(x, ...)
```

## Arguments

- x:

  The rows at which to predict: a numeric matrix, its column names those
  `f` sees, no missing values.

- f:

  The prediction function: a function of a numeric matrix returning one
  finite number per row, or the names of learners of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`"glm"`, `"lasso"`, `"ridge"`, `"enet"`, `"forest"`).

- y:

  The outcome at the rows of `x`: one finite number per row.

- nrep:

  The shuffles of each column: one whole number from 1.

- loss:

  `"se"`, the squared error (the Brier score for a 0/1 outcome), or
  `"ae"`, the absolute error.

- ratio:

  `FALSE` for the increase of the loss, the loss after shuffling less
  the loss before; `TRUE` for their ratio.

- seed:

  The seed of the shuffles: one number, at least 0.

- lctl:

  The settings of the learners named in `f`: a list named by arguments
  of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`fam`, `nfold`, `fold`, `disc`, `alpha`, `nlam`, `nfin`, `ntree`,
  `mtry`, `seed`); empty for a function.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `xpimp`, a list: `imp`, a data frame with one row per
column of `x` (`term`, `imp` the mean increase or ratio, `sd` its
standard deviation over the shuffles); `base`, the loss before
shuffling; `rep`, the increase or ratio of every shuffle (one row per
column, one column per shuffle); `perm`, the shuffles (an array of
`nrow(x)` rows, `nrep` columns and one slice per column of `x`:
`perm[i, r, j]` is the row whose value of column `j` row `i` takes in
shuffle `r`); `loss` and `ratio`.

## Details

`f` is the prediction function: an R function that takes a numeric
matrix with the columns of `x` and returns one prediction per row, or
the names of learners of the package, fitted to `x` and `y` by
[`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
with the settings in `lctl` (one learner is that learner, as
[`rfor()`](https://rootcoder007.github.io/crimlab/reference/rfor.md) or
least squares; several are their stack). The function is called once, on
the rows of `x` followed by every shuffled copy, so a learner is fitted
once. With learners the loss is measured on the rows they were fitted
to; for the loss on new rows, give `f` as a function and `x` and `y` as
the new rows.

The shuffles are drawn in the compiled core from `seed`, the same on
every platform; `perm` returns them, so that every value can be checked.

## References

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
[doi:10.1023/A:1010933404324](https://doi.org/10.1023/A%3A1010933404324)

Fisher, A., Rudin, C. and Dominici, F. (2019). All models are wrong, but
many are useful: learning a variable's importance by studying an entire
class of prediction models simultaneously. *Journal of Machine Learning
Research*, 20(177), 1-81.

## See also

[`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md),
[`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md),
[`xpcp()`](https://rootcoder007.github.io/crimlab/reference/xpcp.md) and
[`xpshap()`](https://rootcoder007.github.io/crimlab/reference/xpshap.md)
for the shape of an input's effect;
[`fair()`](https://rootcoder007.github.io/crimlab/reference/fair.md) and
[`calaud()`](https://rootcoder007.github.io/crimlab/reference/calaud.md)
for the audit of a risk tool's scores.

## Examples

``` r
set.seed(1)
x <- cbind(prior = rpois(200, 2), age = runif(200, 18, 60),
           noise = rnorm(200))
y <- 0.3 * x[, "prior"] - 0.02 * x[, "age"] + rnorm(200, sd = 0.5)
f <- function(m) 0.3 * m[, "prior"] - 0.02 * m[, "age"]
xpimp(x, f, y, nrep = 5)
#> xpimp: permutation importance, 5 shuffles per column; loss se (increase), 0.2852 before shuffling
#>   term       imp         sd
#>  prior 0.3260732 0.02932846
#>    age 0.1047989 0.01306672
#>  noise 0.0000000 0.00000000
xpimp(x, "glm", y, nrep = 5, lctl = list(nfold = 2))
#> xpimp: permutation importance, 5 shuffles per column; loss se (increase), 0.284 before shuffling
#>   term          imp          sd
#>  prior 3.469121e-01 0.031564597
#>    age 1.071655e-01 0.013107077
#>  noise 6.448013e-06 0.000732987
```
