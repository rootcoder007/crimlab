# Posterior probability from a likelihood ratio

Turns a likelihood ratio into the probability of the first hypothesis,
which needs a prior. A large ratio with a small prior can still leave
that probability low; reading the ratio as the probability is the
prosecutor's fallacy. See `proof("lrpost")`.

## Usage

``` r
lrpost(lr, prior)
```

## Arguments

- lr:

  Likelihood ratios: positive numbers.

- prior:

  The prior probability of the first hypothesis, strictly between 0 and
  1: one value, or one per ratio.

## Value

The posterior probabilities.

## Examples

``` r
lrpost(lr = 1e6, prior = 1 / 5e6)
#> [1] 0.1666667
```
