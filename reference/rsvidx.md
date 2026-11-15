# Uniform sample of a stream in one pass

Vitter's Algorithm R: draws `k` items from a stream, each item equally
likely to be retained, in a single pass and using memory proportional to
`k` rather than to the stream's length. `rsvidx()` returns the retained
positions, so the same sample can be applied to a file read line by
line; `rsvsmp()` applies it to a vector in memory.

## Usage

``` r
rsvidx(n, k, seed = 42L)

rsvsmp(x, k, seed = 42L)
```

## Arguments

- n:

  Length of the stream, one non-negative number.

- k:

  Sample size, one non-negative number. Values above the stream's length
  return the whole stream.

- seed:

  Seed for the generator, a whole number from 0 to below \\2^{53}\\.

- x:

  Vector to sample from.

## Value

`rsvidx()` a sorted numeric vector of positions (from 1); `rsvsmp()` the
elements of `x` at them.

## Details

The first `k` items fill the reservoir; item \\i \> k\\ then replaces a
uniformly chosen slot with probability \\k / i\\, so after \\n\\ items
every item is in the sample with probability \\k / n\\ (Vitter 1985).
The positions are returned sorted.

The stream is sampled with the package's own generator seeded by `seed`,
not R's, so the sample is reproducible and R's random stream is left
alone. The generator is C++ `std::mt19937_64`, whose output sequence the
C++ standard fixes, and each draw is Lemire's unbiased
multiply-and-reject on its 64-bit outputs (Lemire 2019), so the same
seed gives the same sample on every platform and C++ library.

Errors: `n` or `k` not one non-negative whole number; `n` or `k` above
2147483647; `seed` not a whole number from 0 to below \\2^{53}\\.

## References

Vitter, J. S. (1985). Random sampling with a reservoir. ACM Transactions
on Mathematical Software 11(1), 37-57.
[doi:10.1145/3147.3165](https://doi.org/10.1145/3147.3165)

Lemire, D. (2019). Fast random integer generation in an interval. ACM
Transactions on Modeling and Computer Simulation 29(1), 3.
[doi:10.1145/3230636](https://doi.org/10.1145/3230636)

## Examples

``` r
rsvidx(1000, 5, seed = 1)
#> [1]  33  67 304 630 638
rsvsmp(letters, 4, seed = 2)
#> [1] "m" "r" "s" "z"
rsvidx(3, 10)
#> [1] 1 2 3
# sample rows of a data frame
df <- data.frame(id = 1:20, v = (1:20)^2)
df[rsvidx(nrow(df), 3, seed = 7), ]
#>   id  v
#> 3  3  9
#> 8  8 64
#> 9  9 81
try(rsvidx(10, -1))
#> Error : rsvidx: `k` must be a single non-negative whole number.
```
