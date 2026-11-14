# Mixture likelihood ratio

The weight of a DNA mixture at one locus for one account of who
contributed against another. Each account names the alleles of the
mixture its known contributors do not explain and the number of unknown
contributors who must explain them. Peak heights, drop-out and
coancestry are not modelled. See `proof("mixlr")`.

## Usage

``` r
mixlr(p, np, xp, nd, xd)
```

## Arguments

- p:

  The population frequencies of the alleles in the mixture: named
  numbers above 0 with a sum of at most 1.

- np, nd:

  The alleles left unexplained under the first and the second account:
  names of `p` (or `character(0)`), at most 20.

- xp, xd:

  The number of unknown contributors under each account.

## Value

A named numeric vector: `num`, `den` and `lr`.

## Details

With no alleles to explain and one unknown, the denominator is the
chance a random person is not excluded from the mixture.

## References

Weir, B. S., Triggs, C. M., Starling, L., Stowell, L. I., Walsh, K. A.
J. and Buckleton, J. (1997). Interpreting DNA mixtures. *Journal of
Forensic Sciences*, 42(2), 213-222.
[doi:10.1520/JFS14100J](https://doi.org/10.1520/JFS14100J)

## Examples

``` r
# victim explains a and b; suspect would explain c and d
p <- c(a = 0.1, b = 0.2, c = 0.05, d = 0.15)
mixlr(p, np = character(0), xp = 0, nd = c("c", "d"), xd = 1)
#>      num      den       lr 
#>  1.00000  0.01500 66.66667 
```
