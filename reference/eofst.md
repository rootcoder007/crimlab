# Empirical orthogonal functions of a space-time matrix

`eofst()` decomposes a matrix of counts or rates, one row per period and
one column per place, into empirical orthogonal functions (EOFs):
spatial patterns (the columns of `eofs`) whose weights over time (the
principal component series, the columns of `pcs`) reproduce the centred
matrix, ordered by the share of its variance each explains. A leading
pattern that loads on a cluster of places, with a series that rises, is
a hot spot that grows; a second pattern of opposite signs is a shift of
crime from one set of places to another.

## Usage

``` r
eofst(z, k = NULL)

# S3 method for class 'eofst'
print(x, ...)
```

## Arguments

- z:

  A numeric matrix (or a data frame of numeric columns) with no missing
  values: one row per period (at least two), one column per place.

- k:

  `NULL` (every component) or the most components to keep: one whole
  number from 1.

- x:

  For [`print()`](https://rdrr.io/r/base/print.html), an `eofst` result.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `eofst`, a list: `eofs` (places by components), `pcs`
(periods by components), `explained` (the share of the variance of each
component), `sdev` and `center` (the column means taken out).

## Details

Each column is centred on its mean over the periods and the centred
matrix `A` is decomposed as `A = U D V'` (Lorenz 1956; Hannachi,
Jolliffe and Stephenson 2007). The EOFs are the columns of `V`, the
series the columns of `A V = U D`, and `explained` is `d^2 / sum(d^2)`;
`sdev` is `d / sqrt(n - 1)`, the standard deviation of each series, as
[`stats::prcomp()`](https://rdrr.io/r/stats/prcomp.html) gives it. An
EOF is defined up to its sign: each is signed so that its largest
loading in absolute value is positive. The compiled kernel finds `V` and
`d` from the eigendecomposition (cyclic Jacobi) of the smaller of `A'A`
and `A A'`; components with `d^2` below `1e-10` of the largest are
dropped, so a matrix of `n` periods has at most `n - 1` components.

## References

Lorenz, E. N. (1956). *Empirical Orthogonal Functions and Statistical
Weather Prediction*. Statistical Forecasting Project, Scientific
Report 1. Massachusetts Institute of Technology.

Hannachi, A., Jolliffe, I. T. and Stephenson, D. B. (2007). Empirical
orthogonal functions and related techniques in atmospheric science: a
review. *International Journal of Climatology*, 27(9), 1119-1152.
[doi:10.1002/joc.1499](https://doi.org/10.1002/joc.1499)

## Examples

``` r
set.seed(3)
yr <- 1:10
z <- outer(yr, c(1, 1, 1, 0, 0, 0)) + matrix(rnorm(60), 10)
colnames(z) <- paste0("ward", 1:6)
e <- eofst(z)
e
#> eofst: 6 components of 10 periods by 6 places
#>  component   sdev explained cumulative
#>          1 5.8950  0.893300     0.8933
#>          2 1.4340  0.052830     0.9461
#>          3 1.0370  0.027640     0.9737
#>          4 0.7186  0.013270     0.9870
#>          5 0.6166  0.009773     0.9968
#>          6 0.3538  0.003218     1.0000
round(e$eofs[, 1], 2)
#> ward1 ward2 ward3 ward4 ward5 ward6 
#>  0.58  0.57  0.57 -0.05 -0.05  0.03 
```
