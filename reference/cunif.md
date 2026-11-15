# Reproducible uniforms (splitmix64)

Returns `n` pseudo-random numbers on \[0, 1) from a 64-bit seed with the
splitmix64 generator. The stream depends only on `seed`, and any program
that implements splitmix64 (Python, Go, a spreadsheet macro) computes
the same numbers, so a subsample or a jitter drawn in one language can
be reproduced exactly in another.

## Usage

``` r
cunif(n, seed)
```

## Arguments

- n:

  The number of uniforms, one whole number from 0.

- seed:

  The seed, one whole number from 0, below \\2^{53}\\ (so that it is
  exact as a double).

## Value

A numeric vector of `n` values in \[0, 1).

## Details

splitmix64 (Steele, Lea and Flood, 2014) advances a 64-bit state by the
constant `0x9E3779B97F4A7C15` and mixes it with two multiply-xorshift
rounds; the top 53 bits of each output divided by \\2^{53}\\ give a
double in \[0, 1). R's own generator and
[`set.seed()`](https://rdrr.io/r/base/Random.html) are neither used nor
changed. It is not a statistical replacement for R's generators.

## References

Steele GL, Lea D, Flood CH (2014). Fast splittable pseudorandom number
generators. *OOPSLA 2014*, 453–472.
[doi:10.1145/2660193.2660195](https://doi.org/10.1145/2660193.2660195)

## Examples

``` r
cunif(3, 42)
#> [1] 0.7415649 0.1599104 0.2786011
all.equal(cunif(3, 7), cunif(10, 7)[1:3])
#> [1] TRUE
set.seed(1); a <- runif(1)
set.seed(1); invisible(cunif(100, 9)); b <- runif(1)
identical(a, b)
#> [1] TRUE
try(cunif(-1, 1))
#> Error : cunif: `n` must be one whole number from 0.
try(cunif(2, 2^53))
#> Error : cunif: `seed` must be one whole number from 0, below 2^53.
```
