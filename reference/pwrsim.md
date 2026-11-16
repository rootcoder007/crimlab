# Power by simulation

The power of a two-arm trial found by simulating it: for each control
arm size in `n` and each effect in `eff`, `nsim` trials are drawn, each
is tested, and the power is the share rejected at `alpha`, with its
Monte Carlo standard error and Wilson's interval (Morris, White and
Crowther 2019). It answers what a formula cannot: an outcome that is
binary or a count, unequal arms, a permutation test.

## Usage

``` r
pwrsim(
  n,
  eff,
  fam = c("normal", "binary", "poisson"),
  sd = 1,
  p0 = 0.5,
  mu0 = 1,
  ratio = 1,
  test = c("welch", "perm"),
  alt = c("two.sided", "greater", "less"),
  alpha = 0.05,
  nsim = 1000,
  nperm = 199,
  level = 0.95,
  seed = NULL
)

# S3 method for class 'pwrsim'
print(x, ...)
```

## Arguments

- n:

  The control arm sizes to try: whole numbers from 2 to 1000000.

- eff:

  The effects to try: differences in means (in the outcome's units; risk
  differences for `"binary"`).

- fam:

  The outcome: `"normal"`, `"binary"` or `"poisson"`.

- sd:

  The standard deviation of a normal outcome: positive.

- p0:

  The control arm's chance of a binary outcome, from 0 to 1; `p0 + eff`
  must be too.

- mu0:

  The control arm's mean of a Poisson outcome: positive; `mu0 + eff`
  must be at least 0.

- ratio:

  The treated arm's size over the control arm's: positive.

- test:

  `"welch"` or `"perm"`.

- alt:

  The alternative: `"two.sided"`, `"greater"` (the treated higher) or
  `"less"`.

- alpha:

  The level of the test: one number strictly between 0 and 1.

- nsim:

  The trials drawn per design: a whole number from 1 to 1000000.

- nperm:

  The relabellings of the permutation test: a whole number from 1 to
  100000 with `1 / (nperm + 1)` at most `alpha`.

- level:

  The level of Wilson's interval: one number strictly between 0 and 1.

- seed:

  `NULL` (drawn from R's generator, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) fixes it) or a
  whole number from 0 to 2147483647.

- x:

  An object of class `pwrsim`.

- ...:

  Not used by [`print()`](https://rdrr.io/r/base/print.html).

## Value

An object of class `pwrsim`, a list: `tab`, a data frame with one row
per design (`n` varying fastest): `n`, `n1` (the treated arm), `eff`,
`power`, `se`, `lo`, `hi` and `rej` (the trials rejected); and `fam`,
`test`, `alt`, `alpha`, `nsim` and `seed`.

## Details

The control arm has `n` units and the treated arm `round(ratio * n)`.
The outcome is normal with standard deviation `sd` and mean 0 in the
control arm and `eff` in the treated arm (`fam = "normal"`), 0 or 1 with
chance `p0` and `p0 + eff` (`"binary"`), or Poisson with mean `mu0` and
`mu0 + eff` (`"poisson"`), drawn by inversion from the package's
generator. The test is Welch's t test (`test = "welch"`; with no spread
in either arm the p-value is 0 when the means differ in the direction
tested, else 1) or the permutation test of the difference in means with
`nperm` relabellings (`"perm"`; `p = (1 + k) / (1 + nperm)`). Each
design has a stream of its own made from `seed` and its place in the
grid, so a design's power does not depend on the others.

## References

Morris, T. P., White, I. R. and Crowther, M. J. (2019). Using simulation
studies to evaluate statistical methods. *Statistics in Medicine*,
38(11), 2074-2102.
[doi:10.1002/sim.8086](https://doi.org/10.1002/sim.8086)

Welch, B. L. (1947). The generalization of 'Student's' problem when
several different population variances are involved. *Biometrika*,
34(1-2), 28-35.
[doi:10.1093/biomet/34.1-2.28](https://doi.org/10.1093/biomet/34.1-2.28)

Wilson, E. B. (1927). Probable inference, the law of succession, and
statistical inference. *Journal of the American Statistical
Association*, 22(158), 209-212.
[doi:10.1080/01621459.1927.10502953](https://doi.org/10.1080/01621459.1927.10502953)

## Examples

``` r
pwrsim(c(20, 50), eff = c(0, 0.5), nsim = 400, seed = 1)
#> pwrsim: normal outcome, welch test (two.sided) at alpha 0.05, 400 trials per design
#>   n n1 eff  power          se         lo         hi rej
#>  20 20 0.0 0.0375 0.009499178 0.02285503 0.06094384  15
#>  50 50 0.0 0.0350 0.009188988 0.02096129 0.05788514  14
#>  20 20 0.5 0.3425 0.023727292 0.29769102 0.39030535 137
#>  50 50 0.5 0.6925 0.023072914 0.64562516 0.73571260 277
pwrsim(100, eff = 0.15, fam = "binary", p0 = 0.3, nsim = 400,
       seed = 2)$tab
#>     n  n1  eff power        se       lo        hi rej
#> 1 100 100 0.15 0.585 0.0246361 0.536129 0.6322539 234
```
