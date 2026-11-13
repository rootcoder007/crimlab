# The Chumbley score for two tool-mark profiles

Compares two striated tool-mark profiles by the deterministic Chumbley
score (Chumbley et al. 2010; Hadler and Morris 2018): each profile is
trimmed and detrended by its lowess smooth, the pair of windows of
`wopt` points with the largest correlation fixes the shift, and windows
of `wval` points at that shift are compared with windows at other shifts
by a standardised Wilcoxon rank-sum statistic.

## Usage

``` r
chmb(a, b, wopt, wval, coarse = 0.25)
```

## Arguments

- a, b:

  The two profiles: finite numbers (heights along the mark).

- wopt, wval:

  The window lengths of the optimisation and validation steps: whole
  numbers from 2 with `wval <= wopt`.

- coarse:

  The lowess span of the detrending: positive.

## Value

A named numeric vector: `nsame`, `ndiff` (window pairs used), `u` (the
statistic) and `p` (its upper normal tail; `NaN` when either set is
empty).

## References

Chumbley, L. S., Morris, M. D., Kreiser, M. J. et al. (2010). Validation
of tool mark comparisons obtained using a quantitative, comparative,
statistical algorithm. *Journal of Forensic Sciences*, 55(4), 953-961.
[doi:10.1111/j.1556-4029.2010.01424.x](https://doi.org/10.1111/j.1556-4029.2010.01424.x)

Hadler, J. R. and Morris, M. D. (2018). An improved version of a tool
mark comparison algorithm. *Journal of Forensic Sciences*, 63(3),
849-855.
[doi:10.1111/1556-4029.13640](https://doi.org/10.1111/1556-4029.13640)

## Examples

``` r
set.seed(1)
m <- cumsum(rnorm(700))
a <- m[1:600] + rnorm(600, sd = 0.2)
b <- m[51:650] + rnorm(600, sd = 0.2)
chmb(a, b, wopt = 120, wval = 30)
#>        nsame        ndiff            u            p 
#> 13.000000000  4.000000000  2.943920289  0.001620417 
```
