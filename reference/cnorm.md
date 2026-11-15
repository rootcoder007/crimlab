# Normal density and log-density (C backend)

`cnormpdf()` evaluates the density of the normal distribution with a
given mean and standard deviation at each value of `x`; `cnormlpdf()`
its logarithm, computed directly from the formula rather than as
`log(dnorm(x))`. Summing log-densities is how a normal model is scored,
and computing them directly keeps the sum finite far into the tails,
where the density itself underflows to zero.

## Usage

``` r
cnormpdf(x, mean = 0, sd = 1)

cnormlpdf(x, mean = 0, sd = 1)
```

## Arguments

- x:

  Numeric vector of quantiles.

- mean:

  The mean \\\mu\\, one finite number.

- sd:

  The standard deviation \\\sigma\\, one positive finite number.

## Value

A numeric vector the length of `x`.

## Details

\$\$f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x -
\mu)^2}{2\sigma^2}\right), \qquad \log f(x) = -\log \sigma - \tfrac12
\log(2\pi) - \frac{(x - \mu)^2}{2\sigma^2}.\$\$ Both agree with
[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html) to rounding.
They are vectorised over `x` only: `mean` must be one finite number and
`sd` one positive finite number. `NA` in `x` gives `NA`; text, a factor
or a list is refused.

## See also

[`stats::dnorm()`](https://rdrr.io/r/stats/Normal.html).

## Examples

``` r
cnormpdf(c(-1, 0, 1))
#> [1] 0.2419707 0.3989423 0.2419707
cnormpdf(c(0, 5, 10), mean = 5, sd = 2)
#> [1] 0.00876415 0.19947114 0.00876415
all.equal(cnormpdf(-2:2), stats::dnorm(-2:2))
#> [1] TRUE
stats::dnorm(50)
#> [1] 0
cnormlpdf(50)
#> [1] -1250.919
y <- c(4.8, 5.1, 5.6, 4.9, 5.3)
sum(cnormlpdf(y, mean = mean(y), sd = stats::sd(y)))
#> [1] -0.9121269
try(cnormpdf(0, sd = 0))
#> Error : cnormpdf: `sd` must be one positive finite number.
try(cnormlpdf(1, mean = c(0, 1)))
#> Error : cnormlpdf: `mean` must be one finite number.
```
