# Accumulated local effects of one input

The effect of one input on a prediction, measured only where the data
are (Apley and Zhu 2020): the range of the column is cut at its
quantiles into intervals, each row's prediction is taken at the upper
and at the lower edge of its own interval with its other inputs as they
are, the differences are averaged within each interval and accumulated
from the lowest edge. Unlike partial dependence
([`xppdp()`](https://rootcoder007.github.io/crimlab/reference/xppdp.md))
no row is moved far from where it lies, so an input correlated with
others (prior arrests and age, say) is not credited with their effects.

## Usage

``` r
xpale(x, f, v, y = NULL, nbin = 10, lctl = list())

# S3 method for class 'xpale'
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

- nbin:

  The intervals asked for: one whole number from 1 to 67108864.

- lctl:

  The settings of the learners named in `f`: a list named by arguments
  of
  [`suplrn()`](https://rootcoder007.github.io/crimlab/reference/suplrn.md)
  (`fam`, `nfold`, `fold`, `disc`, `alpha`, `nlam`, `nfin`, `ntree`,
  `mtry`, `seed`); empty for a function.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `xpale`, a list: `ale`, a data frame with the edges
(`z`) and the centred accumulated effect there (`ale`); `bins`, a data
frame with one row per interval (`lo`, `hi` its edges, `n` its rows,
`eff` the mean difference of the predictions across it); and `term`.

## Details

The edges are the column's minimum and its quantiles at
`1/nbin, 2/nbin, ..., 1`
([`quantile()`](https://rdrr.io/r/stats/quantile.html) type 1), repeated
values dropped, so there are at most `nbin` intervals; each interval
holds the rows above its lower edge up to its upper one, the first its
lower edge too. The effects are centred so that their mean over the
rows, the effect taken as linear within each interval, is 0: the
centring of the authors' 'ALEPlot'. The rows are written and the effects
accumulated in the compiled core; `f` is called once on all the rows.

## References

Apley, D. W. and Zhu, J. (2020). Visualizing the effects of predictor
variables in black box supervised learning models. *Journal of the Royal
Statistical Society Series B*, 82(4), 1059-1086.
[doi:10.1111/rssb.12377](https://doi.org/10.1111/rssb.12377)

## Examples

``` r
set.seed(4)
age <- runif(200, 18, 60)
prior <- rpois(200, (60 - age) / 10)
x <- cbind(prior, age)
f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
xpale(x, f, "age", nbin = 5)
#> xpale: accumulated local effects of age over 5 interval(s)
#>         z         ale
#>  18.12276  0.11761138
#>  26.93682  0.06485629
#>  35.43904  0.01050556
#>  43.48715 -0.03293827
#>  52.64466 -0.06308333
#>  59.85774 -0.07813159
```
