# Bootstrap replicate means (C backend)

`B` resamples of `x`, drawn with replacement and each the length of `x`,
with the mean of each. The spread of these replicate means estimates the
standard error of the sample mean, and their quantiles give a percentile
interval, without assuming the data are normal (Efron and Tibshirani,
1993). The draws come from the core's own 64-bit Mersenne Twister seeded
by `seed`, not from R's generator, so a seed gives the same replicates
in every front end and on every platform, and R's random stream is left
untouched.

## Usage

``` r
cbootmean(x, B = 1000L, seed = 42L)
```

## Arguments

- x:

  Numeric (or logical) vector to resample.

- B:

  The number of replicates, one whole number from 1 to 67108864 (2^26).

- seed:

  The seed of the core's generator, one whole number from 0, below
  \\2^{53}\\.

## Value

A numeric vector of `B` replicate means.

## Details

Replicate \\b\\ draws \\n\\ indices uniformly from \\1, \dots, n\\
(`std::mt19937_64` seeded with `seed`, an unbiased draw by Lemire's
multiply-and-reject) and returns the mean of the values drawn,
accumulated in extended precision. The replicates are one stream: the
first `B` for a seed are the same whatever larger `B` is asked for. An
empty `x` gives `B` values of `NaN`; an `NA` in `x` makes every
replicate that draws it `NA`. A long run can be interrupted.

## References

Efron B, Tibshirani RJ (1993). *An Introduction to the Bootstrap*.
Chapman and Hall.

## Examples

``` r
set.seed(1)
x <- stats::rnorm(50, mean = 5)
reps <- cbootmean(x, B = 500, seed = 7)
c(sample = mean(x), bootstrap = mean(reps))
#>    sample bootstrap 
#>  5.100448  5.097066 
stats::quantile(reps, c(0.025, 0.975))
#>     2.5%    97.5% 
#> 4.864298 5.311439 
identical(cbootmean(x, 10, seed = 1), cbootmean(x, 100, seed = 1)[1:10])
#> [1] TRUE
try(cbootmean(x, B = 0))
#> Error : cbootmean: `B` must be one whole number from 1.
```
