# Column report for a data frame

One row per column: type, missingness, and the summaries that suit the
column's type – mean and standard deviation alongside their robust
counterparts (median, MAD) for a numeric column, and the number of
distinct levels plus the most common one for a categorical column.

## Usage

``` r
profcols(data, qnts = c(0.25, 0.5, 0.75))

# S3 method for class 'profcols'
print(x, ...)
```

## Arguments

- data:

  A data frame whose columns are profiled.

- qnts:

  Numeric vector of quantile probabilities to include for numeric
  columns, one column each named `q` and the percentage (`q25`, `q50`,
  ...); default the quartiles. `NULL` or an empty vector adds none.

- x:

  A result of `profcols()`.

- ...:

  Not used.

## Value

A data frame of class `c("profcols", "data.frame")`, one row per column,
with the columns `column`, `type` (the first class), `n`, `n_missing`,
`pct_missing`, `n_distinct`, `nzero`, `nneg`, `ninf`, `mean`, `sd`,
`median`, `mad`, `min`, `max`, `skewness`, `n_outliers`, `hist` (an
[`inlhist`](https://rootcoder007.github.io/crimlab/reference/inlhist.md)
sketch), `top` and one per quantile; a statistic that does not apply to
the column's type is `NA`. The print method shows the histogram last.

## Details

Reporting the classical and robust centres side by side is the point:
where they disagree, the column has outliers or a heavy tail, and the
mean is not describing it. That is visible in one glance here and in no
single number.

A column is numeric when
[`is.numeric()`](https://rdrr.io/r/base/numeric.html) is `TRUE`; every
other column (character, factor, logical, `Date`, ...) is treated as
categorical, and `top` is its most frequent non-missing value as text
(ties broken by the order of `sort(table())`). For a numeric column,
over its non-missing values: `sd` is the sample standard deviation
(\\n - 1\\), `mad` the median absolute deviation scaled by 1.4826,
`skewness` the population moment ratio \\m_3 / m_2^{3/2}\\, `n_outliers`
the count outside the Tukey fences (1.5 interquartile ranges beyond the
quartiles), and the quantiles are type 7. Infinite values are counted in
`ninf` (and in `nneg`, `n_distinct`) and left out of every other numeric
statistic, which are taken over the finite values only: a column with an
infinite value keeps a finite mean and range, and a column with no
finite value has `NA` statistics. Errors when `data` is not a data
frame, has no columns, or has duplicate or empty column names.

## See also

[`inlhist`](https://rootcoder007.github.io/crimlab/reference/inlhist.md),
[`freqtbl`](https://rootcoder007.github.io/crimlab/reference/freqtbl.md).

## Examples

``` r
set.seed(9)
df <- data.frame( # simulated columns
  clean = stats::rnorm(100),
  skewed = c(stats::rnorm(99), 500),
  grade = sample(c("a", "b", "c"), 100, TRUE),
  gappy = c(rep(NA, 10), stats::runif(90))
)

profcols(df)
#> Column profile: 4 columns
#>  column      type   n n_missing pct_missing n_distinct nzero nneg ninf     mean
#>   clean   numeric 100         0           0        100     0   54    0 -0.05351
#>  skewed   numeric 100         0           0        100     0   48    0  4.87406
#>   grade character 100         0           0          3    NA   NA   NA       NA
#>   gappy   numeric 100        10          10         90     0    0    0  0.50195
#>       sd   median    mad       min      max skewness n_outliers  top     q25
#>   0.9588 -0.17079 0.8777 -2.617706   2.6820  0.31449          3 <NA> -0.7457
#>  50.0224  0.04449 0.9517 -2.737955 500.0000  9.84352          1 <NA> -0.7947
#>       NA       NA     NA        NA       NA       NA         NA    c      NA
#>   0.2917  0.52428 0.3714  0.007639   0.9953 -0.04396          0 <NA>  0.2545
#>       q50    q75       hist
#>  -0.17079 0.4220 ▁▁▃█▇█▃▃▂▁
#>   0.04449 0.6215 █        ▁
#>        NA     NA       <NA>
#>   0.52428 0.7339 █▆▆█▆▇█▆█▆

# The mean and the median agree on `clean` and disagree sharply on
# `skewed`, which is the outlier announcing itself.
p <- profcols(df)
p[p$column %in% c("clean", "skewed"), c("column", "mean", "median")]
#> Column profile: 2 columns
#>  column     mean   median
#>   clean -0.05351 -0.17079
#>  skewed  4.87406  0.04449

# The histogram column shows shape no summary number carries.
p[, c("column", "hist")]
#> Column profile: 4 columns
#>  column       hist
#>   clean ▁▁▃█▇█▃▃▂▁
#>  skewed █        ▁
#>   grade       <NA>
#>   gappy █▆▆█▆▇█▆█▆

# Zeros, negatives and infinities are counted separately, because
# each breaks a different downstream computation (a log, a square
# root, an average).
p[, c("column", "nzero", "nneg", "ninf")]
#> Column profile: 4 columns
#>  column nzero nneg ninf
#>   clean     0   54    0
#>  skewed     0   48    0
#>   grade    NA   NA   NA
#>   gappy     0    0    0

# deciles instead of quartiles
names(profcols(df, qnts = c(0.1, 0.9)))
#>  [1] "column"      "type"        "n"           "n_missing"   "pct_missing"
#>  [6] "n_distinct"  "nzero"       "nneg"        "ninf"        "mean"       
#> [11] "sd"          "median"      "mad"         "min"         "max"        
#> [16] "skewness"    "n_outliers"  "hist"        "top"         "q10"        
#> [21] "q90"        

# a table with no columns has nothing to profile
try(profcols(df[, 0]))
#> Error : profcols: `data` has no columns to profile.
```
