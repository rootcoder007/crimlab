# Multivariate Hawkes process: types that excite each other

Fits a multivariate Hawkes process with exponential kernels to event
times of several types (offence types, gangs, areas): \$\$\lambda_c(t) =
\mu_c + \beta \sum_l A\_{cl} \sum\_{t_j \< t, \\ type_j = l} e^{-\beta
(t - t_j)},\$\$ where \\A\_{cl}\\ is the expected number of events of
type `c` triggered by one event of type `l` (retaliation between gangs,
one crime type leading to another). The rates, the matrix and the shared
decay `beta` (or a given one) are fitted by maximum likelihood with the
analytic gradient (Ozaki's recursion, O(n k)) by BFGS over their
logarithms. The process is stationary when the spectral radius `rad` of
\\A\\ is below 1; with one type it is the exponential process of
[`hwkfit()`](https://rootcoder007.github.io/crimlab/reference/hwkfit.md)
with `eta` \\= A\\.

## Usage

``` r
mvhwk(times, type, hor, beta = NULL, start = NULL, maxit = 500)
```

## Arguments

- times:

  Event times: numbers sorted in `[0, hor]`, at least 10.

- type:

  The type of each event: labels, at most 10 types.

- hor:

  The end of the window: one positive number.

- beta:

  `NULL` to estimate the decay, or a given one: one positive number.

- start:

  `NULL` or the starting parameters in the order of `par` (positive).

- maxit:

  The most BFGS steps: one whole number from 1.

## Value

An object of class `mvhwk`, a list: `par`, `se` and `cov` (named
`mu.<type>`, `a.<to>.<from>`, `beta`); `mu`, the rates; `a`, the matrix
(rows the triggered type, columns the triggering type); `beta`; `rad`;
`loglik`, `aic`; `conv`, `iter`; `n` and `hor`.

## References

Hawkes, A. G. (1971). Spectra of some self-exciting and mutually
exciting point processes. *Biometrika*, 58(1), 83-90.
[doi:10.1093/biomet/58.1.83](https://doi.org/10.1093/biomet/58.1.83)

Linderman, S. W. and Adams, R. P. (2014). Discovering latent network
structure in point process data. *Proceedings of the 31st International
Conference on Machine Learning*, 32(2), 1413-1421.

## Examples

``` r
a <- hwksim(c(log(0.3), 0.3, 1), 400, seed = 1)
b <- hwksim(c(log(0.2), 0.4, 1), 400, seed = 2)
d <- rbind(data.frame(t = a$t, type = "a"), data.frame(t = b$t, type = "b"))
d <- d[order(d$t), ]
mvhwk(d$t, d$type, 400)
#> mvhwk: 302 events of 2 types; beta 0.9417787, spectral radius 0.2974369
#> background rates:
#>         a         b 
#> 0.3202407 0.2003798 
#> excitation (row: triggered, column: triggering):
#>              a          b
#> a 2.794820e-01 0.06558192
#> b 5.088176e-09 0.29743689
```
