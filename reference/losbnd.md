# Mean time in custody from completed and current stays

The mean length of stay of the people released (or of the cases
disposed) leaves out everyone still inside. Their elapsed time is a
lower bound on their stay, so the mean over everyone admitted is at
least `lower`, the mean of the completed stays and the elapsed times
together; the released-only mean falls short of it by `m / (n + m)`
times the excess of the mean elapsed time over the released mean, and
understates the truth for certain when that excess is positive. No upper
bound exists without a cap on how long the current stays will last. See
`proof("losbnd")`.

## Usage

``` r
losbnd(done, open)
```

## Arguments

- done:

  Completed stays: non-negative numbers.

- open:

  Elapsed time of the stays still open: non-negative numbers.

## Value

A named numeric vector: `mean` (of the completed stays), `age` (the mean
elapsed time of the open ones), `lower` (the lower bound on the mean
over all) and `bias` (`lower - mean`, the least by which the completed
mean understates when it is positive).

## References

Manski, C. F. (2003). *Partial Identification of Probability
Distributions*. Springer.
[doi:10.1007/b97478](https://doi.org/10.1007/b97478)

## Examples

``` r
losbnd(done = c(30, 45, 60, 90, 120), open = c(100, 150, 200))
#>    mean     age   lower    bias 
#>  69.000 150.000  99.375  30.375 
```
