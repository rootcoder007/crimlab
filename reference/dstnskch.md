# Distinct-value count in fixed memory

HyperLogLog: estimates how many distinct values a column holds using a
fixed register per bucket, regardless of the column's length or
cardinality. `dstnskch()` builds or extends the registers, `dstncnt()`
turns them into an estimate, and `skchmrg()` combines two sketches.

## Usage

``` r
dstnskch(x, p = 14L, regstr = NULL)

dstncnt(regstr)

skchmrg(a, b)
```

## Arguments

- x:

  A vector, compared as text.

- p:

  Log2 of the register count, a whole number from 4 to 20 (default 14).

- regstr:

  Registers from a previous call (with the same `p`), to extend; `NULL`
  starts afresh. For `dstncnt()`, the registers.

- a, b:

  Register sets to merge.

## Value

`dstnskch()` and `skchmrg()` an integer vector of \\2^p\\ registers;
`dstncnt()` the estimated number of distinct values.

## Details

Each value is hashed and the first `p` bits pick one of \\m = 2^p\\
registers; the register keeps the longest run of leading zeros (plus
one) seen in the remaining bits. The estimate is the bias-corrected
harmonic mean \\\alpha_m m^2 / \sum_j 2^{-M_j}\\ of Flajolet et al.
(2007).

The estimate carries a relative standard error of about
`1.04 / sqrt(2^p)`, so 0.8% at the default. It is an estimate:
`length(unique(x))` is exact when the column fits in memory. Below
roughly `2.5 * 2^p` distinct values the estimator switches to linear
counting, which is near-exact in that range.

Values are hashed with SHA-256 (FIPS 180-4) of their UTF-8 text, so the
sketch is identical on every platform and across sessions. Values are
compared as text
([`as.character()`](https://rdrr.io/r/base/character.html)), so `1`,
`1.0` and `"1"` are one value; `NA` is skipped. Since a register only
keeps a maximum, folding chunks with `regstr` and merging separate
sketches with `skchmrg()` give the same registers as one sketch of all
the data.

Errors: `p` not a whole number from 4 to 20; `regstr` not of length
\\2^p\\; registers that are not whole numbers from 0 to 64; two sketches
of different sizes.

## References

Flajolet, P., Fusy, E., Gandouet, O. and Meunier, F. (2007).
HyperLogLog: the analysis of a near-optimal cardinality estimation
algorithm. Analysis of Algorithms 2007, 137-156.

## Examples

``` r
set.seed(1)
x <- sample(1:5000, 200000, replace = TRUE) # simulated
dstncnt(dstnskch(x))
#> [1] 5039.269
length(unique(x))
#> [1] 5000
dstncnt(dstnskch(c("a", "b", "c", "a", "b")))
#> [1] 3.000275
a <- dstnskch(x[1:100000])
b <- dstnskch(x[100001:200000])
dstncnt(skchmrg(a, b))
#> [1] 5039.269
try(skchmrg(dstnskch(1:10, p = 4), dstnskch(1:10, p = 5)))
#> Error : skchmrg: both sketches must have the same number of registers.
```
