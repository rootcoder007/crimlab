# Capture-recapture estimates of a hidden population

How many offenders, victims or incidents no list records, from the
overlap of two or three lists (police records, hospital records, a
survey). With two lists the Lincoln-Petersen estimate `n1 n2 / m`
assumes the lists are independent; when being on one list changes the
odds of being on the other by at most a factor `kappa`, the truth lies
within `[n1 n2 / (kappa m), kappa n1 n2 / m]` (`lower`, `upper`), and no
data can narrow that box. Chapman's estimate \`(n1 + 1)(n2 + 1) / (m +
1)

- 1\` is nearly unbiased, with Seber's standard error.

## Usage

``` r
crcap(x, kappa = 1, int = NULL)
```

## Arguments

- x:

  `c(n1, n2, m)`, the sizes of two lists and their overlap (named so or
  in that order), or the seven cells of three lists, named by
  membership: `"100"`, `"010"`, `"001"`, `"110"`, `"101"`, `"011"`,
  `"111"`. Counts from 0.

- kappa:

  The largest factor by which being on one list changes the odds of
  being on the other: one number from 1 (two lists).

- int:

  `NULL`, or three-way log interactions at which to give the missing
  cell (three lists).

## Value

An object of class `crcap`, a list: `type` (`"two"` or `"three"`);
`est`, a named vector (two lists: `petersen`, `lower`, `upper`,
`chapman`, `se`, `floor`; three: `observed`, `missing`, `total`);
`pairs` (three lists: `lists`, `n1`, `n2`, `m`, `floor`, `petersen`,
`chapman`) and `sens` (with `int`: `int`, `missing`, `total`), otherwise
`NULL`.

## Details

With three lists the seven observed cells identify the never-listed cell
only under an assumption on the three-way interaction of a log-linear
model: with none, it is `m100 m010 m001 m111 / (m110 m101 m011)`, and a
three-way log interaction `t` multiplies it by `exp(-t)`; `int` gives
that sensitivity. The pairwise estimates are given beside it, each at
least the number seen on either of its two lists (`floor`). See
`proof("crcap")`.

## References

Chapman, D. G. (1951). Some properties of the hypergeometric
distribution with applications to zoological sample censuses.
*University of California Publications in Statistics*, 1, 131-160.

Seber, G. A. F. (1970). The effects of trap response on tag recapture
estimates. *Biometrics*, 26(1), 13-22.
[doi:10.2307/2529040](https://doi.org/10.2307/2529040)

Fienberg, S. E. (1972). The multiple recapture census for closed
populations and incomplete 2^k contingency tables. *Biometrika*, 59(3),
591-603.
[doi:10.1093/biomet/59.3.591](https://doi.org/10.1093/biomet/59.3.591)

## Examples

``` r
crcap(c(n1 = 500, n2 = 400, m = 80), kappa = 1.5)
#> crcap: two lists 
#>  petersen     lower     upper   chapman        se     floor 
#> 2500.0000 1666.6667 3750.0000 2479.2593  224.0262  820.0000 
crcap(c("100" = 120, "010" = 90, "001" = 60, "110" = 30, "101" = 20,
        "011" = 15, "111" = 8), int = c(-0.5, 0, 0.5))
#> crcap: three lists 
#> observed  missing    total 
#>      343      576      919 
#>  lists  n1  n2  m floor petersen  chapman
#>    1-2 178 143 38   283 669.8421 659.9231
#>    1-3 178 103 28   253 654.7857 640.9310
#>    2-3 143 103 23   223 640.3913 623.0000
#>   int  missing     total
#>  -0.5 949.6635 1292.6635
#>   0.0 576.0000  919.0000
#>   0.5 349.3617  692.3617
```
