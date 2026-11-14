# Risk terrain model

Risk terrain modelling (Caplan, Kennedy and Miller 2011) relates the
crimes counted in the cells of a grid to risk layers, the presence or
density of features (bars, transit stops, vacant lots) around each cell,
by Poisson regression. The exponentiated coefficient of a layer is its
relative risk value; the cells' relative risk combines the layers.

## Usage

``` r
rtmfit(y, lay, cell = NULL, rad = NULL, kind = "prox", level = 0.95)
```

## Arguments

- y:

  The crime counts of the cells: numbers from 0, not all 0.

- lay:

  The layers: columns with one value per cell (`cell` `NULL`), or a list
  of feature coordinates, a two-column matrix for each layer. Names of
  the columns or list elements name the layers.

- cell:

  `NULL`, or the cells' centres: a matrix of two columns with one row
  per cell.

- rad:

  With `cell`: the distance within which features count, or the kernel
  bandwidth: one positive number.

- kind:

  With `cell`: `"prox"`, `"dens"` or `"kde"`.

- level:

  The confidence level of the intervals: one number strictly between 0
  and 1.

## Value

A list: `coef`, a data frame with one row per term (`term`, `est` the
coefficient, `se`, `z`, `p`, `rrv` the relative risk value `exp(est)`
and its interval `lo`, `hi`); `cells`, a data frame with one row per
cell (`rr` the relative risk, the product of the layers' relative risk
values to the power of the cell's layer values; `score`, `rr` over its
smallest value; `fitted`; `y`); `lay`, the matrix of layers; `dev` the
deviance; `aic`; `conv`, whether the fit converged; `level`.

## Details

The layers are given either as columns (`cell` `NULL`: a matrix, a data
frame or a list with one value per cell each) or built here from the
coordinates of the features: `cell` holds the cells' centres and `lay` a
list of feature coordinates (one two-column matrix per layer), and a
cell's value is, by `kind`, 1 when a feature lies within `rad` and 0
otherwise (`"prox"`), the number of features within `rad` (`"dens"`), or
the Gaussian kernel density of the features with bandwidth `rad`,
standardised over the cells to mean 0 and standard deviation 1 (`"kde"`;
0 everywhere when it does not vary). The fit is maximum likelihood by
iteratively reweighted least squares, as
[`stats::glm()`](https://rdrr.io/r/stats/glm.html) with
`family = poisson()` does; a layer that the others explain is aliased,
its coefficient `NA`.

## References

Caplan, J. M., Kennedy, L. W. and Miller, J. (2011). Risk terrain
modeling: brokering criminological theory and GIS methods for crime
forecasting. *Justice Quarterly*, 28(2), 360-381.
[doi:10.1080/07418825.2010.486037](https://doi.org/10.1080/07418825.2010.486037)

## Examples

``` r
rtmfit(c(1, 2, 4, 8, 3, 5), list(bars = c(0, 1, 0, 1, 0, 1),
                                 stops = c(0, 0, 1, 1, 1, 0)))
#> $coef
#>          term       est        se         z          p      rrv        lo
#> 1 (Intercept) 0.2876821 0.5188745 0.5544348 0.57928134 1.333333 0.4822497
#> 2        bars 0.9162907 0.4557327 2.0105880 0.04436899 2.500000 1.0233392
#> 3       stops 0.9162907 0.4557327 2.0105880 0.04436899 2.500000 1.0233392
#>         hi
#> 1 3.686426
#> 2 6.107457
#> 3 6.107457
#> 
#> $cells
#>     rr score   fitted y
#> 1 1.00  1.00 1.333333 1
#> 2 2.50  2.50 3.333333 2
#> 3 2.50  2.50 3.333333 4
#> 4 6.25  6.25 8.333333 8
#> 5 2.50  2.50 3.333333 3
#> 6 2.50  2.50 3.333333 5
#> 
#> $lay
#>      bars stops
#> [1,]    0     0
#> [2,]    1     0
#> [3,]    0     1
#> [4,]    1     1
#> [5,]    0     1
#> [6,]    1     0
#> 
#> $dev
#> [1] 1.609242
#> 
#> $aic
#> [1] 25.89929
#> 
#> $conv
#> [1] TRUE
#> 
#> $level
#> [1] 0.95
#> 
# layers from feature locations
cell <- as.matrix(expand.grid(x = 1:4, y = 1:4))
bars <- rbind(c(1, 1), c(3, 4))
y <- c(3, 1, 0, 0, 2, 1, 0, 1, 0, 0, 1, 2, 0, 1, 2, 3)
rtmfit(y, list(bars = bars), cell = cell, rad = 1.5, kind = "dens")$coef
#>          term       est       se         z         p       rrv         lo
#> 1 (Intercept) -1.791759 1.000000 -1.791760 0.0731715 0.1666667 0.02347725
#> 2        bars  2.261763 1.030776  2.194233 0.0282187 9.6000000 1.27313009
#>          hi
#> 1  1.183179
#> 2 72.388518
```
