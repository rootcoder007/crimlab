# Two-level model for multivariate trace evidence

`tlvc()` estimates the background of the two-level normal model of
Aitken and Lucy (2004) from a reference collection,
`x_ij = theta_i + e_ij`, items `theta_i ~ N(mu, C)` and measurements
within an item `e_ij ~ N(0, U)`: the overall mean, the pooled
within-item covariance and the ANOVA between-item covariance (unbiased
for any pattern of replication). `tlvlr()` gives the likelihood ratio of
same source against different sources for two items (a control and a
recovered fragment, say) from their mean vectors, computed in logarithms
from Cholesky factors so that atypical or distant items give a finite
value.

## Usage

``` r
tlvc(x, item)

tlvlr(xa, xb, na, nb, bg)
```

## Arguments

- x:

  The reference measurements: a numeric matrix, one row per observation
  and one column per variable.

- item:

  The item of each row: labels, every item at least once, some item more
  than once.

- xa, xb:

  The mean vectors of the two items compared: one number per variable.

- na, nb:

  Their numbers of replicates: positive.

- bg:

  The background, as `tlvc()` returns it.

## Value

`tlvc()`: a list of `mu` (the mean), `u` and `c` (the within- and
between-item covariance matrices; `c` need not be positive definite,
which `tlvlr()` then refuses). `tlvlr()`: a named vector, `ln` and
`log10`, the natural and base-10 logarithms of the ratio.

## References

Aitken, C. G. G. and Lucy, D. (2004). Evaluation of trace evidence in
the form of multivariate data. *Journal of the Royal Statistical
Society, Series C*, 53(1), 109-122.
[doi:10.1046/j.0035-9254.2003.05271.x](https://doi.org/10.1046/j.0035-9254.2003.05271.x)

## Examples

``` r
set.seed(1)
th <- matrix(rnorm(20, 10, 2), 10, 2)
x <- th[rep(1:10, each = 3), ] + matrix(rnorm(60, 0, 0.2), 30, 2)
bg <- tlvc(x, rep(1:10, each = 3))
tlvlr(th[1, ], th[1, ] + 0.05, 3, 3, bg)
#>       ln    log10 
#> 5.978629 2.596485 
tlvlr(th[1, ], th[2, ], 3, 3, bg)
#>        ln     log10 
#> -387.3127 -168.2078 
```
