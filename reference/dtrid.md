# Which deterrence effects a design identifies

For the response `y = b0 + b1 p + b2 s + b3 c` to the certainty `p`,
severity `s` and celerity (delay) `c` of the sanction regimes observed,
a dimension that never varies has no identified coefficient, and the
coefficients are identified exactly when the design `(1, p, s, c)` has
full rank. Gives the rank and, for each dimension, whether its effect is
identified given the others (full rank, or the rank falls without it).
With `log = TRUE` the design is `(1, log p, log s, c)`: the expected
cost of an offence `p s exp(-rho c)` is linear on that scale, and with a
response `y` the ratio of the celerity to the certainty coefficient,
`-b3 / b1`, is the discount rate `rho` the response implies.

## Usage

``` r
dtrid(p, s, c, y = NULL, log = FALSE, rho = NULL)
```

## Arguments

- p, s, c:

  The certainty, severity and delay of each regime: finite numbers, one
  per regime (`p` and `s` positive with `log = TRUE`).

- y:

  `NULL`, or the response of each regime: finite numbers (the design
  must then have full rank).

- log:

  `TRUE` for the log design, `FALSE` for the plain one.

- rho:

  `NULL`, or a discount rate at which to give the expected cost of each
  regime: one finite number.

## Value

A list: `rank`; `ident` and `const`, for `certainty`, `severity` and
`celerity`, whether the effect is identified and whether the dimension
is constant; `coef`, `NULL` or the least-squares coefficients; `disc`,
`NULL`, or the implied discount rate (`NA` but for the log design);
`cost`, `NULL` or `p * s * exp(-rho * c)`.

## References

Nagin, D. S. (2013). Deterrence in the twenty-first century. *Crime and
Justice*, 42(1), 199-263.
[doi:10.1086/670398](https://doi.org/10.1086/670398)

## Examples

``` r
# certainty varied at fixed severity and celerity
dtrid(p = c(0.1, 0.3, 0.5), s = c(2, 2, 2), c = c(30, 30, 30))
#> $rank
#> [1] 2
#> 
#> $ident
#> certainty  severity  celerity 
#>      TRUE     FALSE     FALSE 
#> 
#> $const
#> certainty  severity  celerity 
#>     FALSE      TRUE      TRUE 
#> 
#> $coef
#> NULL
#> 
#> $disc
#> NULL
#> 
#> $cost
#> NULL
#> 
set.seed(2)
p <- runif(40, 0.1, 0.9)
s <- runif(40, 1, 10)
c <- runif(40, 0, 12)
y <- 2 - 0.8 * log(p) - 0.8 * log(s) + 0.08 * c + rnorm(40, sd = 0.1)
dtrid(p, s, c, y, log = TRUE)$disc
#> [1] 0.09283505
```
