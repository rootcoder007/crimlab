# Full pairwise correlation table

Computes the correlation of every pair of numeric columns and returns it
in long form – one row per pair, which is easier to sort, filter and
join than a matrix – with the number of complete pairs each figure rests
on.

## Usage

``` r
cortbl(data, mthd = c("spearman", "pearson"), minprs = 3L)

# S3 method for class 'cortbl'
print(x, ...)
```

## Arguments

- data:

  A data frame; its numeric columns are correlated in pairs.

- mthd:

  `"spearman"` (default) or `"pearson"`.

- minprs:

  Minimum complete pairs required before a correlation is computed, a
  whole number of at least 2 (default 3). Below it the pair is `NA`
  rather than a number computed from almost nothing.

- x:

  A result of `cortbl()`.

- ...:

  Not used.

## Value

A data frame of class `c("cortbl", "data.frame")` with one row per pair
of numeric columns and the columns `x` and `y` (character, the column
names), `correlation` (numeric, in \\\[-1, 1\]\\ or `NA`) and `nprs`
(integer); the method is in `attr(, "method")`.

## Details

Non-numeric columns are ignored. For each pair, rows where either value
is missing are dropped (pairwise deletion), so different pairs can rest
on different rows; `nprs` says how many. Spearman's rank correlation
(the default) is the Pearson correlation of the ranks, robust to
outliers and to monotone transformations; see
[`topcors`](https://rootcoder007.github.io/crimlab/reference/topcors.md)
for why it is the default. A pair with fewer than `minprs` complete
rows, or with a constant column (no variance), gets `NA`. The pairs are
listed in column order
([`utils::combn()`](https://rdrr.io/r/utils/combn.html)).

Errors: `data` not a data frame; fewer than two numeric columns;
`minprs` not a whole number of at least 2; an invalid `mthd`.

## See also

[`topcors`](https://rootcoder007.github.io/crimlab/reference/topcors.md)
for just the strongest.

## Examples

``` r
set.seed(1)
df <- data.frame(a = stats::rnorm(100), b = stats::rnorm(100))
df$c <- df$a + stats::rnorm(100, sd = 0.2) # simulated

cortbl(df)
#> Correlations (spearman): 3 pairs
#>  x y correlation nprs
#>  a b      0.0387  100
#>  a c      0.9665  100
#>  b c      0.0151  100

# Pearson instead of Spearman
cortbl(df, mthd = "pearson")
#> Correlations (pearson): 3 pairs
#>  x y correlation nprs
#>  a b   -0.000994  100
#>  a c    0.974707  100
#>  b c   -0.012038  100

# nprs shows how much data each figure rests on.
gappy <- df
gappy$a[1:80] <- NA
cortbl(gappy)
#> Correlations (spearman): 3 pairs
#>  x y correlation nprs
#>  a b      0.1654   20
#>  a c      0.9594   20
#>  b c      0.0151  100

# A pair with too little overlap is NA, not a number from nothing.
thin <- df
thin$a[1:99] <- NA
cortbl(thin)$correlation
#> [1]         NA         NA 0.01508551

# one numeric column cannot be correlated
try(cortbl(data.frame(a = 1:5, b = letters[1:5])))
#> Error : cortbl: need at least two numeric columns to correlate.
```
