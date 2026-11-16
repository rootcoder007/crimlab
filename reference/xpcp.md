# Ceteris paribus profiles of single rows

The prediction for one person as one input moves over a grid while the
person's other inputs stay as they are: what the risk tool would have
said had only that input been different. One profile per row of `newx`;
over all the rows of `x` the profiles are the individual conditional
expectation curves of Goldstein et al. (2015), and their mean is the
partial dependence of
[`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md).
A steep profile on a protected or proxy input is a decision that turns
on that input alone.

## Usage

``` r
xpcp(
  x,
  f,
  v,
  newx = x[1L, , drop = FALSE],
  y = NULL,
  grid = 20,
  gtype = c("even", "quant"),
  lctl = list()
)

# S3 method for class 'xpcp'
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

- v:

  The input: the name or the number of one column of `x`.

- newx:

  The rows to profile: a numeric matrix with the columns of `x` (by
  default its first row).

- y:

  The outcome at the rows of `x`, needed only to fit the learners named
  in `f`.

- grid:

  The number of grid points, one whole number from 2, or the values of
  the grid (at least two).

- gtype:

  `"even"` for points evenly spaced over the column's range, `"quant"`
  for its quantiles.

- lctl:

  The settings of the learners named in `f`: a list named by arguments
  of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`fam`, `nfold`, `fold`, `disc`, `alpha`, `nlam`, `nfin`, `ntree`,
  `mtry`, `seed`); empty for a function.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `xpcp`, a list: `prof`, the predictions (one row per
row of `newx`, one column per grid value); `grid`; `pred`, the
prediction at each row of `newx` as it is; and `term`.

## Details

The grid is taken from the column of `x`, as in
[`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md).
The rows are written in the compiled core; `f` is called once on all of
them.

## References

Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015). Peeking
inside the black box: visualizing statistical learning with plots of
individual conditional expectation. *Journal of Computational and
Graphical Statistics*, 24(1), 44-65.
[doi:10.1080/10618600.2014.907095](https://doi.org/10.1080/10618600.2014.907095)

Biecek, P. and Burzykowski, T. (2021). *Explanatory Model Analysis*.
Chapman and Hall/CRC, New York, ch. 10.

## Examples

``` r
set.seed(3)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpcp(x, f, "prior", newx = x[1:2, ], grid = c(0, 2, 4, 6))
#> xpcp: ceteris paribus profiles of 2 row(s) over 4 values of prior
#>       pred         lo        hi
#>  0.1186177 0.07546754 0.6211458
#>  0.2891560 0.08321185 0.6457738
```
