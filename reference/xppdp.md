# Partial dependence of a prediction function on one input

The average prediction when one input is set to a value for every row
and the others keep their values (Friedman 2001, sec. 8.2): at each
value `g` of a grid, the mean of the predictions with column `v` set to
`g` in every row. It is the effect of the input averaged over the data,
as a risk tool's user would read it; where the input is correlated with
the others it averages over combinations the data never hold, which
[`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md)
avoids.

## Usage

``` r
xppdp(x, f, v, y = NULL, grid = 20, gtype = c("even", "quant"), lctl = list())

# S3 method for class 'xppdp'
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

An object of class `xppdp`, a list: `pd`, a data frame with the grid
(`grid`) and the partial dependence there (`pd`); `term`, the input's
name; and `n`, the rows averaged.

## Details

The grid is `grid` values of the column, evenly spaced from its minimum
to its maximum (`gtype = "even"`) or at its quantiles
(`gtype = "quant"`,
[`quantile()`](https://rdrr.io/r/stats/quantile.html)'s default type 7,
repeated values dropped), or the values given. The rows are written in
the compiled core; `f` is called once on all of them.

## References

Friedman, J. H. (2001). Greedy function approximation: a gradient
boosting machine. *Annals of Statistics*, 29(5), 1189-1232.
[doi:10.1214/aos/1013203451](https://doi.org/10.1214/aos/1013203451)

## See also

[`xpcp()`](https://rootcoder007.github.io/crimlab/reference/xpcp.md) for
the same curve for single rows,
[`xpale()`](https://rootcoder007.github.io/crimlab/reference/xpale.md).

## Examples

``` r
set.seed(2)
x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xppdp(x, f, "age", grid = 5)
#> xppdp: partial dependence on age, averaged over 100 rows
#>      grid        pd
#>  18.57225 0.3701719
#>  28.72549 0.3092730
#>  38.87872 0.2541271
#>  49.03196 0.2055543
#>  59.18520 0.1638529
```
