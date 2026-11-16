# Sampling Shapley values of single predictions

Shares out one person's prediction among the inputs: the Shapley value
of input `j` is its contribution averaged over every order in which the
inputs could be revealed, the inputs not yet revealed taking values from
the rows of `x` (the reference data). Strumbelj and Kononenko (2014)
estimate it by sampling: each sample draws an order of the inputs and a
row `w` of `x`, and takes the prediction at the person's values on `j`
and the inputs before it (`w`'s on the rest) less the prediction with
`j` too taken from `w`. The value is the mean over `nsamp` samples, its
standard error their standard deviation over `sqrt(nsamp)`. Over all the
inputs the values add up, in expectation, to the person's prediction
less the mean prediction over `x`.

## Usage

``` r
xpshap(
  x,
  f,
  newx = x[1L, , drop = FALSE],
  y = NULL,
  nsamp = 200,
  seed = 1,
  lctl = list()
)

# S3 method for class 'xpshap'
print(x, ...)
```

## Arguments

- x:

  The reference data: a numeric matrix, its column names those `f` sees,
  no missing values (the rows the learners in `f` are fitted to).

- f:

  The prediction function: a function of a numeric matrix returning one
  finite number per row, or the names of learners of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`"glm"`, `"lasso"`, `"ridge"`, `"enet"`, `"forest"`).

- newx:

  The rows to explain: a numeric matrix with the columns of `x` (by
  default its first row).

- y:

  The outcome at the rows of `x`, needed only to fit the learners named
  in `f`.

- nsamp:

  The samples of each value: one whole number from 1.

- seed:

  The seed of the draws: one number, at least 0.

- lctl:

  The settings of the learners named in `f`: a list named by arguments
  of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`fam`, `nfold`, `fold`, `disc`, `alpha`, `nlam`, `nfin`, `ntree`,
  `mtry`, `seed`); empty for a function.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `xpshap`, a list: `phi` and `se`, the values and
their standard errors (one row per row of `newx`, one column per input);
`pred`, the prediction at each row of `newx`; `base`, the mean
prediction over the rows of `x`; and `nsamp`.

## Details

The orders and rows are drawn in the compiled core from `seed`, the same
on every platform; `f` is called once on all the rows.

## References

Strumbelj, E. and Kononenko, I. (2014). Explaining prediction models and
individual predictions with feature contributions. *Knowledge and
Information Systems*, 41(3), 647-665.
[doi:10.1007/s10115-013-0679-x](https://doi.org/10.1007/s10115-013-0679-x)

Shapley, L. S. (1953). A value for n-person games. In Kuhn, H. W. and
Tucker, A. W. (eds), *Contributions to the Theory of Games II*, 307-317.
Princeton University Press.

## Examples

``` r
set.seed(5)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpshap(x, f, newx = x[1:2, ], nsamp = 100)
#> xpshap: sampling Shapley values, 100 samples each; mean prediction 0.2695
#>        pred       prior         age
#>  0.09888651 -0.10661348 -0.06044138
#>  0.27926231  0.07142368 -0.05804795
```
