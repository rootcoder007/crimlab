# Describing a table and keeping its categories whole

Profiles of a data frame, frequency and correlation tables, and guards
against categories that change silently.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Describing a table

### Column report for a data frame

[`profcols()`](https://rootcoder007.github.io/crimlab/reference/profcols.md):
One row per column: type, missingness, and the summaries that suit the
column’s type - mean and standard deviation alongside their robust
counterparts (median, MAD) for a numeric column, and the number of
distinct levels plus the most common one for a categorical column.

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
#>  column      type   n n_missing pct_missing n_distinct nzero nneg ninf     mean      sd
#>   clean   numeric 100         0           0        100     0   54    0 -0.05351  0.9588
#>  skewed   numeric 100         0           0        100     0   48    0  4.87406 50.0224
#>   grade character 100         0           0          3    NA   NA   NA       NA      NA
#>   gappy   numeric 100        10          10         90     0    0    0  0.50195  0.2917
#>    median    mad       min      max skewness n_outliers  top     q25      q50    q75
#>  -0.17079 0.8777 -2.617706   2.6820  0.31449          3 <NA> -0.7457 -0.17079 0.4220
#>   0.04449 0.9517 -2.737955 500.0000  9.84352          1 <NA> -0.7947  0.04449 0.6215
#>        NA     NA        NA       NA       NA         NA    c      NA       NA     NA
#>   0.52428 0.3714  0.007639   0.9953 -0.04396          0 <NA>  0.2545  0.52428 0.7339
#>        hist
#>  ▁▁▃█▇█▃▃▂▁
#>  █        ▁
#>        <NA>
#>  █▆▆█▆▇█▆█▆

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
#>  [1] "column"      "type"        "n"           "n_missing"   "pct_missing" "n_distinct" 
#>  [7] "nzero"       "nneg"        "ninf"        "mean"        "sd"          "median"     
#> [13] "mad"         "min"         "max"         "skewness"    "n_outliers"  "hist"       
#> [19] "top"         "q10"         "q90"

# a table with no columns has nothing to profile
try(profcols(df[, 0]))
#> Error : profcols: `data` has no columns to profile.
```

### Inline histogram for a numeric vector

[`inlhist()`](https://rootcoder007.github.io/crimlab/reference/inlhist.md):
Draws a one-line sketch of a column’s distribution, as block characters
(or ASCII where the console cannot render them). The point is density of
information: a mean and a standard deviation cannot tell you a column is
bimodal, skewed or has a gap, and this can, in the width of a table
cell. `profcols` shows one for every numeric column.

``` r

set.seed(1)
# A symmetric distribution peaks in the middle (simulated values).
inlhist(stats::rnorm(1000))
#> [1] "▁▂▄██▇▄▂▁▁"

# A skewed one leans left.
inlhist(stats::rexp(1000))
#> [1] "█▄▂▁▁▁▁▁▁▁"

# Bimodality is visible here and in no single summary number.
inlhist(c(stats::rnorm(500, -3), stats::rnorm(500, 3)))
#> [1] "▁▅█▄▁▁▃█▆▁"

# A constant column has no spread to show.
inlhist(rep(5, 10))
#> [1] "     █    "

# more bins for more detail
inlhist(stats::rexp(1000), bins = 30)
#> [1] "█▇▅▄▃▃▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁ ▁ ▁  ▁ ▁"

# text is not a distribution
try(inlhist(c("a", "b")))
#> Error : inlhist: `x` must be numeric, got character.
```

### Frequency table for one column

[`freqtbl()`](https://rootcoder007.github.io/crimlab/reference/freqtbl.md):
Counts how often each value of a column occurs, with percentages of all
rows and of the non-missing rows, and missing values counted as their
own row. It is the first look at a categorical variable (a charge type,
a region, a response category).

``` r

df <- data.frame(grade = c("a", "b", "b", "c", NA, "b"))
freqtbl(df, "grade")
#> Frequency table: 6 values
#>  value n      pct pctvalid                  bar
#>      b 3 50.00000       60 ####################
#>      a 1 16.66667       20              #######
#>      c 1 16.66667       20              #######
#>   <NA> 1 16.66667       NA              #######

# The two percentage columns differ exactly by the missingness.
freqtbl(df, "grade")[, c("pct", "pctvalid")]
#>       pct pctvalid
#>  50.00000       60
#>  16.66667       20
#>  16.66667       20
#>  16.66667       NA

# A vector works directly.
freqtbl(c(1, 1, 2, 3, 3, 3))
#> Frequency table: 6 values
#>  value n      pct pctvalid                  bar
#>      3 3 50.00000 50.00000 ####################
#>      1 2 33.33333 33.33333        #############
#>      2 1 16.66667 16.66667              #######

# Sorted value order rather than frequency order.
freqtbl(c("c", "a", "b", "a"), sort = FALSE)
#> Frequency table: 4 values
#>  value n pct pctvalid                  bar
#>      a 2  50       50 ####################
#>      b 1  25       25           ##########
#>      c 1  25       25           ##########

# Long tails are folded so the percentages still total 100.
set.seed(1)
ft <- freqtbl(sample(letters, 500, TRUE), maxlvl = 5)
ft
#> Frequency table: 500 values
#>    value   n  pct pctvalid                  bar
#>        y  25  5.0      5.0                    #
#>        e  24  4.8      4.8                    #
#>        t  24  4.8      4.8                    #
#>        a  23  4.6      4.6                    #
#>        f  23  4.6      4.6                    #
#>  (other) 381 76.2     76.2 ####################
sum(ft$pct)
#> [1] 100

# a data frame needs a column name
try(freqtbl(df))
#> Error : freqtbl: `col` is required when `data` is a data frame.
```

### Full pairwise correlation table

[`cortbl()`](https://rootcoder007.github.io/crimlab/reference/cortbl.md):
Computes the correlation of every pair of numeric columns and returns it
in long form - one row per pair, which is easier to sort, filter and
join than a matrix - with the number of complete pairs each figure rests
on.

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

### Strongest pairwise correlations in a data frame

[`topcors()`](https://rootcoder007.github.io/crimlab/reference/topcors.md):
Ranks the numeric column pairs by the strength of their association and
returns the strongest, so a wide table’s structure can be read without
squinting at a correlation matrix. Use it on a new data set to find
near-duplicate columns, derived columns and the relationships worth a
closer look.

``` r

set.seed(1)
n <- 200
df <- data.frame(
  a = stats::rnorm(n),
  b = stats::rnorm(n),
  grade = sample(letters[1:3], n, TRUE)
)
df$c <- df$a * 2 + stats::rnorm(n, sd = 0.1) # strongly related to a
df$d <- exp(df$a) # monotone but curved

topcors(df)
#> Top correlations (spearman)
#>  x y correlation abscor
#>  a d      1.0000 1.0000
#>  a c      0.9984 0.9984
#>  c d      0.9984 0.9984
#>  a b     -0.0616 0.0616
#>  b d     -0.0616 0.0616
#>  b c     -0.0597 0.0597

# The curved pair is ranked correctly by Spearman and understated by
# Pearson, which is why Spearman is the default.
tc <- topcors(df, mthd = "spearman")
tc[tc$x == "a" & tc$y == "d", "correlation"]
#> [1] 1
tp <- topcors(df, mthd = "pearson")
tp[tp$x == "a" & tp$y == "d", "correlation"]
#> [1] 0.8496382

# Filter to the pairs worth looking at, and keep the top two.
topcors(df, n = 2, minabs = 0.5)
#> Top correlations (spearman)
#>  x y correlation abscor
#>  a d       1.000  1.000
#>  a c       0.998  0.998

# a threshold above 1 is an error
try(topcors(df, minabs = 2))
#> Error : topcors: `minabs` must be a single value in [0, 1].
```

### Drop empty or constant columns and rows

[`drpempty()`](https://rootcoder007.github.io/crimlab/reference/drp.md),
[`drpconst()`](https://rootcoder007.github.io/crimlab/reference/drp.md):
[`drpempty()`](https://rootcoder007.github.io/crimlab/reference/drp.md)
removes rows or columns that are entirely missing.
[`drpconst()`](https://rootcoder007.github.io/crimlab/reference/drp.md)
removes columns that hold a single distinct value. Use them after
reading a spreadsheet export, which often carries blank padding rows,
empty columns and columns that repeat one value (a year, a jurisdiction)
on every row.

``` r

df <- data.frame(
  keep = c(1, 2, NA),
  allna = c(NA, NA, NA),
  constant = c(7, 7, 7)
)

drpempty(df)
#> -- Empty rows and columns dropped
#>    rows 3 -> 3, columns 3 -> 2
#>    columns dropped (all missing): allna
#>   keep constant
#> 1    1        7
#> 2    2        7
#> 3   NA        7
attr(drpempty(df), "dropped")
#> -- Empty rows and columns dropped
#>    rows 3 -> 3, columns 3 -> 2
#>    columns dropped (all missing): allna
#> [1] "allna"

drpconst(df)
#> -- Constant columns dropped
#>    columns 3 -> 1
#>    columns dropped (one value): allna, constant
#>   keep
#> 1    1
#> 2    2
#> 3   NA
attr(drpconst(df), "dropped")
#> -- Constant columns dropped
#>    columns 3 -> 1
#>    columns dropped (one value): allna, constant
#> [1] "allna"    "constant"

# Rows only.
drpempty(data.frame(a = c(1, NA), b = c(2, NA)), which = "rows")
#> -- Empty rows and columns dropped
#>    rows 2 -> 1, columns 2 -> 2
#>   a b
#> 1 1 2

# A column of NAs plus one value is constant by default, and not when
# NA is treated as a value of its own.
x <- data.frame(v = c(NA, NA, 5))
ncol(drpconst(x))
#> -- Constant columns dropped
#>    columns 1 -> 0
#>    columns dropped (one value): v
#> [1] 0
ncol(drpconst(x, naasval = TRUE))
#> -- Constant columns dropped
#>    columns 1 -> 1
#> [1] 1

# only data frames
try(drpempty(c(1, NA)))
#> Error : drpempty: `data` must be a data frame.
```

### Tidy the column names of an ingested table

[`clncolnms()`](https://rootcoder007.github.io/crimlab/reference/clncolnms.md):
Converts column names to a consistent, syntactically valid form:
transliterated to ASCII, non-alphanumerics collapsed to a single
separator, and duplicates disambiguated with a numeric suffix. Use it
right after reading an open-data file, so the rest of the script can
refer to columns without backticks.

``` r

clncolnms(c("Total  Population (2021)", "% change", "Ville / City",
            "dup", "dup"))
#> -- Column names cleaned
#>    Total  Population (2021) -> total_population_2021
#>    % change -> pct_change
#>    Ville / City -> ville_city
#>    dup -> dup_2
#> [1] "total_population_2021" "pct_change"            "ville_city"           
#> [4] "dup"                   "dup_2"

# Applied to a data frame, with the original names retained.
df <- data.frame(`Total Pop` = 1:2, `% change` = 3:4, check.names = FALSE)
cleaned <- clncolnms(df)
#> -- Column names cleaned
#>    Total Pop -> total_pop
#>    % change -> pct_change
names(cleaned)
#> [1] "total_pop"  "pct_change"
attr(cleaned, "orignms")
#> [1] "Total Pop" "% change"

# Other cases.
clncolnms("Total Pop", case = "lower_camel")
#> -- Column names cleaned
#>    Total Pop -> totalPop
#> [1] "totalPop"
clncolnms("Total Pop", case = "upper_camel")
#> -- Column names cleaned
#>    Total Pop -> TotalPop
#> [1] "TotalPop"
clncolnms("Total Pop", case = "screaming_snake")
#> -- Column names cleaned
#>    Total Pop -> TOTAL_POP
#> [1] "TOTAL_POP"

# a different separator, and a name that starts with a digit
clncolnms(c("2021 count", "Rate per 1000"), sep = ".")
#> -- Column names cleaned
#>    2021 count -> x2021.count
#>    Rate per 1000 -> rate.per.1000
#> [1] "x2021.count"   "rate.per.1000"

# an unknown case is an error
try(clncolnms("a b", case = "kebab"))
#> Error : clncolnms: `case` must be one of "snake", "lower_camel", "upper_camel", "screaming_snake", "none", not "kebab".
```

### Multivariate outliers by Mahalanobis distance

[`mahaoutl()`](https://rootcoder007.github.io/crimlab/reference/mahaoutl.md):
Finds the rows of a numeric table that are unusual as a combination of
values, not just on one variable at a time, and gives each row a
distance, a p-value and an outlier flag. Use it before modelling to find
records that are probably data-entry or linkage errors.

``` r

set.seed(1)
n <- 200
df <- data.frame(height = stats::rnorm(n, 170, 10)) # simulated
df$weight <- df$height * 0.5 + stats::rnorm(n, 0, 5)

# A row that is ordinary on each variable but unlikely jointly:
# short and heavy.
df[1, ] <- list(height = 155, weight = 100)

out <- mahaoutl(df)
head(out, 3)
#> Mahalanobis outliers (robust): 1 of 3 rows beyond alpha = 0.001
#>  row distance   p_value logpvalue outlier
#>    1    4.301 9.611e-05    -9.250    TRUE
#>   32    2.900 1.492e-02    -4.205   FALSE
#>   14    2.807 1.944e-02    -3.940   FALSE

# It is flagged, although both values lie inside the other rows' ranges.
out$row[1] == 1
#> [1] TRUE
range(df$height[-1])
#> [1] 147.8530 194.0162
range(df$weight[-1])
#> [1]  67.07546 101.82631

# A cluster of bad rows drags the mean toward itself; the median and
# MAD move less, so the robust distance of such a row is the larger.
bad <- df
bad[2:21, ] <- list(height = 150, weight = 110)
dclassical <- mahaoutl(bad, robust = FALSE)
drobust <- mahaoutl(bad, robust = TRUE)
dclassical$distance[dclassical$row == 2] <
  drobust$distance[drobust$row == 2]
#> [1] TRUE

# A duplicated column has no distance defined, and is refused.
dup <- df
dup$height2 <- dup$height
try(mahaoutl(dup))
#> Error : mahaoutl: the columns are collinear, so no Mahalanobis distance is defined. Drop the redundant column(s) -- drpconst() and topcors() will find them.
```

### Runs of consecutive missing values

[`missruns()`](https://rootcoder007.github.io/crimlab/reference/missruns.md):
Finds the maximal stretches of consecutive `NA` in each column, with
where each begins and how long it is.

``` r

# One long outage and two isolated gaps, with the same total count.
df <- data.frame(
  outage = c(1, 2, NA, NA, NA, NA, 7, 8),
  scattered = c(1, NA, 3, 4, NA, 6, NA, NA)
)
sum(is.na(df$outage)) == sum(is.na(df$scattered))
#> [1] TRUE

missruns(df)
#> Runs of consecutive missing values: 2 
#>     column start end length
#>     outage     3   6      4
#>  scattered     7   8      2

# Isolated gaps too, by lowering the threshold.
missruns(df, minrun = 1)
#> Runs of consecutive missing values: 4 
#>     column start end length
#>     outage     3   6      4
#>  scattered     7   8      2
#>  scattered     2   2      1
#>  scattered     5   5      1

# A complete column has no runs.
missruns(data.frame(x = 1:5))
#> Runs of consecutive missing values: 0

try(missruns(df, minrun = 0))
#> Error : missruns: `minrun` must be a single positive integer.
```

### Text map of where the missing values are

[`missgrid()`](https://rootcoder007.github.io/crimlab/reference/missgrid.md):
Draws the missingness of a whole table as a grid, one character per cell
block, a console map that needs no graphics device, so it works over a
terminal connection, in a log, and in a plain-text summary.

``` r

set.seed(1)
df <- data.frame(
  complete = 1:100,
  block = c(rep(NA, 30), 31:100),
  scattered = ifelse(stats::runif(100) < 0.3, NA, 1), # simulated gaps
  mostlygone = c(1:10, rep(NA, 90))
)
missgrid(df)
#> ── Missingness map ───────────────────────────────────────────────
#>          cbsm
#>          olco
#>          moas
#>          pctt
#>          lktl
#>          e ey
#>          t rg
#>          e eo
#>            dn
#>             e
#>        1  █▒ 
#>        6  █░ 
#>       11  █▒█
#>       16  █ █
#>       21  █▒█
#>       26  █░█
#>       31   ░█
#>       36   ░█
#>       41    █
#>       46   ░█
#>       51   ▒█
#>       56   ░█
#>       61   ░█
#>       66   ▒█
#>       71    █
#>       76    █
#>       81    █
#>       86   ▓█
#>       91   ▒█
#>       96    █
#> 
#>   legend: ' ' none  '█' all missing
#> ──────────────────────────────────────────────────────────────────

# a coarser map with short labels
missgrid(df, height = 5, width = 4)
#> ── Missingness map ───────────────────────────────────────────────
#>          cbsm
#>          olco
#>          moas
#>          pctt
#>        1  █░▒
#>       21  ▒░█
#>       41   ░█
#>       61   ░█
#>       81   ░█
#> 
#>   legend: ' ' none  '█' all missing
#> ──────────────────────────────────────────────────────────────────

try(missgrid(df, height = 0))
#> Error : missgrid: `height` must be a single positive integer.
```

### Which columns are missing together

[`misspat()`](https://rootcoder007.github.io/crimlab/reference/misspat.md):
Counts the distinct PATTERNS of missingness across rows - which
combination of columns is missing in each row - rather than the
per-column rates `profcols` reports. Use it before deciding how to
handle gaps (drop rows, impute, model separately), because the right
choice depends on whether columns are missing together.

``` r

# Two columns missing in the SAME rows: one structural gap.
structural <- data.frame(
  id = 1:10,
  a = c(rep(NA, 3), 4:10),
  b = c(rep(NA, 3), 4:10)
)
misspat(structural)
#> Missingness patterns, columns in pattern order: id, a, b 
#>  pattern n_rows pctrows n_missing columns
#>      ...      7      70         0        
#>      .XX      3      30         2    a, b

# The same per-column rates, but missing independently.
scattered <- data.frame(
  id = 1:10,
  a = c(rep(NA, 3), 4:10),
  b = c(1:7, rep(NA, 3))
)
misspat(scattered)
#> Missingness patterns, columns in pattern order: id, a, b 
#>  pattern n_rows pctrows n_missing columns
#>      ...      4      40         0        
#>      ..X      3      30         1       b
#>      .X.      3      30         1       a

# A complete frame has exactly one pattern.
misspat(data.frame(x = 1:3, y = 4:6))
#> Missingness patterns, columns in pattern order: x, y 
#>  pattern n_rows pctrows n_missing columns
#>       ..      3     100         0

# keep only the most common pattern
misspat(scattered, maxpat = 1)
#> Missingness patterns, columns in pattern order: id, a, b 
#>  pattern n_rows pctrows n_missing columns
#>      ...      4      40         0

# a frame without columns has no patterns
try(misspat(data.frame()))
#> Error : misspat: `data` has no columns.
```

### Missingness in one line per question

[`misssmry()`](https://rootcoder007.github.io/crimlab/reference/misssmry.md):
The scalar summaries of missingness: how much of the table is missing,
how many rows are complete, and how many columns have any or only
missing values. Use it as a first look at a new extract, or to report
data completeness in a methods section; `misspat` then shows which
columns are missing together.

``` r

df <- data.frame(a = c(1, NA, 3), b = c(NA, NA, 3), c = 1:3)
misssmry(df)
#>             n_rows             n_cols          n_missing        pct_missing 
#>            3.00000            3.00000            3.00000           33.33333 
#>    n_complete_rows  pct_complete_rows n_cols_any_missing n_cols_all_missing 
#>            1.00000           33.33333            2.00000            0.00000

# A complete table is all zeros but for its dimensions.
misssmry(data.frame(x = 1:3, y = 4:6))
#>             n_rows             n_cols          n_missing        pct_missing 
#>                  3                  2                  0                  0 
#>    n_complete_rows  pct_complete_rows n_cols_any_missing n_cols_all_missing 
#>                  3                100                  0                  0

# a wholly empty column
ms <- misssmry(data.frame(x = 1:4, note = NA))
ms[c("pct_missing", "n_cols_all_missing")]
#>        pct_missing n_cols_all_missing 
#>                 50                  1

# a vector is not a table
try(misssmry(c(1, NA, 3)))
#> Error : misssmry: `data` must be a data frame.
```

### Little’s test for data missing completely at random

[`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md):
Tests the MCAR assumption: that whether a value is missing is unrelated
to any value in the data, observed or not.

``` r

set.seed(1)
n <- 300
x <- stats::rnorm(n) # simulated
y <- x + stats::rnorm(n)

# Missing completely at random: a coin flip decides, so the test
# should not reject.
mcar <- data.frame(x = x, y = y)
mcar$y[sample(n, 90)] <- NA
mcartst(mcar)
#> Little's MCAR test: no departure from MCAR detected
#>   Statistic 1.094 on 1 df, p = 0.296
#>   Patterns 2, variables 2, rows used 300, EM iterations 13

# Missing depending on the OTHER, observed variable: not MCAR, and
# detectable, because the pattern's mean of x is shifted.
mar <- data.frame(x = x, y = y)
mar$y[x > 0.4] <- NA
mcartst(mar)
#> Little's MCAR test: MCAR rejected: the missingness is related to the data
#>   Statistic 180.118 on 1 df, p = 4.57e-41
#>   Patterns 2, variables 2, rows used 300, EM iterations 69

# Complete data has one pattern and nothing to test.
mcartst(data.frame(a = x, b = y))$df
#> [1] 0

# A duplicated column makes the likelihood degenerate, and is refused.
dup <- mcar
dup$x2 <- dup$x
try(mcartst(dup))
#> Error : mcartst: the columns are collinear (or one is constant), so the multivariate normal likelihood Little's test is built on is degenerate. Drop the redundant column(s) first -- topcors() and drpconst() will find them.
```

### Multiple imputation by predictive mean matching

[`mipmm()`](https://rootcoder007.github.io/crimlab/reference/mipmm.md):
Fills the missing values of numeric columns `m` times by predictive mean
matching, each completed data set from its own chain of chained
equations, and pools an analysis of the completed data sets by Rubin’s
rules: the step after
[`mcartst()`](https://rootcoder007.github.io/crimlab/reference/mcartst.md)
has said how the values are missing. Predictive mean matching imputes
only values that were observed, so an imputed count stays a whole number
and an imputed age stays in its range (van Buuren 2018, sec. 3.4).

``` r

set.seed(11)
n <- 200
d <- data.frame(age = round(runif(n, 18, 60)), male = rbinom(n, 1, 0.5))
d$priors <- rpois(n, exp(-1 + 0.04 * d$age))
d$days <- round(20 + 3 * d$priors + 5 * d$male + rnorm(n, sd = 5))
d$priors[sample(n, 40)] <- NA
d$days[sample(n, 30)] <- NA
r <- mipmm(d, c("priors", "days"), c("age", "male"), fy = "days",
           fx = c("priors", "male"), seed = 1)
r
#> mipmm: 70 missing value(s) in 2 column(s), 5 imputations by predictive mean matching
#>         term    est     se    df     lo     hi         p    riv    fmi
#>  (Intercept) 18.629 0.8058 40.95 17.002 20.257 4.065e-25 0.3591 0.2977
#>       priors  3.432 0.2773 27.11  2.863  4.001 1.154e-12 0.5178 0.3849
#>         male  6.242 0.8145 61.78  4.614  7.871 1.550e-10 0.2468 0.2227
head(r$imp$priors)
#>    [,1] [,2] [,3] [,4] [,5]
#> 11    2    1    0    1    0
#> 17    0    1    2    0    2
#> 18    1    0    0    3    2
#> 19    0    0    0    0    2
#> 21    1    1    3    1    1
#> 23    0    1    0    0    0
```

### Duplicated rows, with their groups

[`duprows()`](https://rootcoder007.github.io/crimlab/reference/duprows.md):
Returns the rows that share a combination of `cols` with at least one
other row, grouped so the duplicates sit together. A duplicate is
usually a join that fanned out or a re-release appended rather than
replaced, both of which are visible only once the offending rows are in
front of you.

``` r

df <- data.frame(
  id = c(1, 2, 2, 3, 3, 3),
  value = c("a", "b", "b", "c", "d", "c")
)

# Duplicated on every column.
duprows(df)
#> -- Duplicated rows
#>    4 row(s) of 6
#>   id value dupe_count
#> 1  2     b          2
#> 2  2     b          2
#> 3  3     c          2
#> 4  3     c          2

# Duplicated on the identifier alone, which catches more.
duprows(df, "id")
#> -- Duplicated rows
#>    5 row(s) of 6
#>   id value dupe_count
#> 1  2     b          2
#> 2  2     b          2
#> 3  3     c          3
#> 4  3     d          3
#> 5  3     c          3

# No duplicates gives zero rows, not an error.
duprows(data.frame(x = 1:3))
#> -- Duplicated rows
#>    0 row(s) of 3
#> [1] x          dupe_count
#> <0 rows> (or 0-length row.names)

# a column that does not exist is an error
try(duprows(df, "code"))
#> Error : duprows: no such column(s): code.
```

### Fit and rank candidate distributions

[`dstfit()`](https://rootcoder007.github.io/crimlab/reference/dstfit.md):
Fits candidate distributions to one numeric variable by maximum
likelihood, ranks them by AIC, and reports for each the BIC, a
Kolmogorov-Smirnov distance, an Anderson-Darling statistic and
percentile bootstrap intervals for its parameters. A Hill estimate of
the tail index is added for positive data. Use it to choose the model of
a count of offences, a sentence length or a time to reoffending before
modelling it.

``` r

set.seed(1)
x <- stats::rgamma(300, shape = 2, rate = 0.5) # simulated
f <- dstfit(x, nboot = 20)
f
#> Distributions fitted by maximum likelihood: 300 values
#>  rank      family npar loglik  aic  bic      ks      ad
#>     1     weibull    2 -673.5 1351 1359 0.02543  0.1965
#>     2       gamma    2 -674.3 1353 1360 0.03906  0.5089
#>     3   lognormal    2 -691.7 1387 1395 0.07364  3.4610
#>     4 exponential    1 -709.4 1421 1424 0.14763 13.0561
#>     5           t    3 -711.2 1428 1440 0.08261  3.4834
#>     6      normal    2 -717.2 1438 1446 0.09560  4.6589
#>     7      pareto    2 -945.7 1895 1903 0.36066 73.1955
#>                                       params
#>                 shape = 1.523, scale = 4.349
#>                shape = 2.0060, rate = 0.5125
#>               meanlog = 1.095, sdlog = 0.812
#>                                rate = 0.2555
#>  location = 3.686, scale = 2.283, df = 8.140
#>                     mean = 3.914, sd = 2.642
#>                xmin = 0.1680, alpha = 0.3474
#> Hill tail index (k = 17): 5.66 (2.97 to 8.35)
f$fits$gamma
#>     shape      rate 
#> 2.0061190 0.5125443

# simulated counts with extra zeros
n <- ifelse(stats::runif(300) < 0.3, 0, stats::rpois(300, 4))
dstfit(n, nboot = 0)$table[, c("family", "aic")]
#>    family      aic
#> 1     zip 1183.895
#> 2    zinb 1185.895
#> 3  negbin 1280.316
#> 4  normal 1382.649
#> 5       t 1384.649
#> 6 poisson 1448.838

# a family the data do not allow
try(dstfit(c(-1, 2, 3, 4, 5), fams = "gamma"))
#> Error : dstfit: `gamma` needs every value above 0.
```

### Do distinct groups differ? A heterogeneity test before pooling

[`dstntst()`](https://rootcoder007.github.io/crimlab/reference/dstntst.md):
Tests whether the estimates of one indicator for distinct groups (First
Nations, Metis and Inuit; Aboriginal and Torres Strait Islander people;
any groups with published intervals) differ more than their intervals
allow, and says how much of their spread is real: the evidence for
reporting each group distinctly rather than as one pooled figure.

``` r

# your own table: made-up numbers, with the pooled row named
own <- data.frame(
  group = c("Aboriginal", "Torres Strait Islander", "Both"),
  est = c(24.0, 31.0, 24.6),
  lo = c(22.1, 26.2, 22.8),
  hi = c(25.9, 35.8, 26.4)
)
dstntst(own, "group", "est", "lo", "hi", pooled = "Both")
#> Distinctions test: 2 groups
#>   Cochran's Q 7.06 on 1 df, p = 0.00787
#>   I-squared: 86% of the variation is beyond chance
#>   Pooled: 24.95 (se 0.901, inverse-variance mean)
#>   Left out: Both
#>                   group estimate        se
#>              Aboriginal       24 0.9694056
#>  Torres Strait Islander       31 2.4490246

# simulated victimisation rates of three police divisions
div <- data.frame(
  division = c("North", "Central", "South"),
  rate = c(4.1, 5.3, 4.4), lo = c(3.6, 4.7, 3.9), hi = c(4.6, 5.9, 4.9)
)
dstntst(div, "division", "rate", "lo", "hi")
#> Distinctions test: 3 groups
#>   Cochran's Q 9.42 on 2 df, p = 0.00899
#>   I-squared: 79% of the variation is beyond chance
#>   Pooled: 4.52 (se 0.155, inverse-variance mean)
#>    group estimate        se
#>    North      4.1 0.2551067
#>  Central      5.3 0.3061281
#>    South      4.4 0.2551067
```

## Categorical integrity

### Recode categories by name

[`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
[`deccodes()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md),
[`relbl()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md):
[`grdrcd()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md)
recodes a categorical column by an explicit mapping from old to new
labels and refuses any observed value the mapping does not name: silent
pass-through is how group labels get corrupted.
[`deccodes()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md)
decodes numeric category codes with a code-to-label dictionary (codes as
names), the step that turns the codes of an imported file into its
categories.
[`relbl()`](https://rootcoder007.github.io/crimlab/reference/grdrcd.md)
relabels a variable by name, never by position, and returns a factor
whose levels follow the mapping.

``` r

# simulated race codes of an imported file
codes <- c(1, 2, 2, 3, 1, NA)
deccodes(codes, c("1" = "White", "2" = "Black", "3" = "Indigenous"))
#> [1] "White"      "Black"      "Black"      "Indigenous" "White"      NA          
#> attr(,"recode")
#> attr(,"recode")$mapping
#>            1            2            3 
#>      "White"      "Black" "Indigenous" 
#> 
#> attr(,"recode")$kept
#> character(0)
grdrcd(c("M", "F", "F", "X"), c(M = "Male", F = "Female"), keep = "X")
#> [1] "Male"   "Female" "Female" "X"     
#> attr(,"recode")
#> attr(,"recode")$mapping
#>        M        F 
#>   "Male" "Female" 
#> 
#> attr(,"recode")$kept
#> [1] "X"
relbl(factor(c("a", "b", "a")), c(a = "Adult", b = "Youth"))
#> [1] Adult Youth Adult
#> attr(,"recode")
#> attr(,"recode")$mapping
#>     a     b 
#> Adult Youth 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
#> Levels: Adult Youth
try(grdrcd(c("M", "F", "U"), c(M = "Male", F = "Female")))
#> Error : grdrcd: values with no mapping: 'U'. Every observed category must be mapped by name (or listed in `keep`).
```

### A factor with declared, verified levels

[`grdlvls()`](https://rootcoder007.github.io/crimlab/reference/grdlvls.md):
Builds a factor whose levels are exactly the ones declared, in the
declared order, and refuses any value outside them. A declared reference
must be the first level, so the group every model contrast is measured
against is stated rather than left to alphabetical order.

``` r

grdlvls(c("Male", "Female", "Male"), c("Male", "Female"), ref = "Male")
#> [1] Male   Female Male  
#> Levels: Male Female
try(grdlvls(c("Male", "Female", "X"), c("Male", "Female")))
#> Error : grdlvls: values outside the declared levels: 'X'.
try(grdlvls("Male", c("Male", "Female"), ref = "Female"))
#> Error : grdlvls: the declared reference 'Female' is not the first level ('Male'); put the reference first in `lvls`.
```

### Decode a labelled import by its value labels

[`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md):
Decodes a column of codes by the value labels it carries (the `labels`
attribute that readers of SPSS, Stata and SAS files attach, a named
vector label = code), or by `vlabs` (code = label) when it carries none.
The codes alone are not the categories: decoding by code keeps every
observation with its own label whatever order the labels are listed in.

``` r

x <- structure(c(2, 1, 3, 2), labels = c(Female = 1, Male = 2, Other = 3))
declbld(x)
#> [1] Male   Female Other  Male  
#> attr(,"recode")
#> attr(,"recode")$mapping
#>      1      2      3 
#> Female   Male  Other 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
#> Levels: Female Male Other
declbld(c(10, 20, 10), c("20" = "Remand", "10" = "Sentenced"))
#> [1] Sentenced Remand    Sentenced
#> attr(,"recode")
#> attr(,"recode")$mapping
#>        20        10 
#>    Remand Sentenced 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
#> Levels: Sentenced Remand
try(declbld(c(1, 2)))
#> Error : declbld: `x` carries no `labels` attribute and `vlabs` was not given; the codes alone are not the categories.
```

### Refuse a categorical column where a 0/1 indicator is required

[`grdbin()`](https://rootcoder007.github.io/crimlab/reference/grdbin.md):
Checks that a column meant as a numeric 0/1 indicator (a treatment, an
outcome) is one. A factor or text column is refused rather than
converted: [`as.numeric()`](https://rdrr.io/r/base/numeric.html) on a
factor returns the level indices (1, 2, …), not the data, and a
mis-ordered level relabels every observation. Encode it explicitly
first, for example `as.integer(x == "treated")`.

``` r

grdbin(c(0, 1, 1, NA), "custody")
try(grdbin(factor(c("no", "yes")), "custody"))
#> Error : grdbin: column 'custody' is categorical ('no', 'yes'...). It is not coerced: as.numeric() on a factor returns level indices, not the data. Encode it first, for example as.integer(x == "treated").
try(grdbin(c(0, 1, 2), "custody"))
#> Error : grdbin: column 'custody' must be 0 or 1 (saw: 0, 1, 2).
```

### Prove a recode with a before/after cross-tabulation

[`vrfrcd()`](https://rootcoder007.github.io/crimlab/reference/vrfrcd.md):
Cross-tabulates the original and the recoded column and checks that the
recode is a function of the original label (no category was split
between two new ones), that it follows the declared mapping (a value the
mapping does not name must be unchanged), that no row was lost or added
and that missingness did not change. This is how swapped groups are
caught before any model runs.

``` r

o <- c("M", "F", "F", NA)
m <- c(M = "Male", F = "Female")
vrfrcd(o, grdrcd(o, m), m)
try(vrfrcd(o, c("Female", "Male", "Male", NA), m))
#> Error : vrfrcd: 'F' was mapped to 'Male' but the declared mapping says 'Female'. THIS is how groups get swapped; fix the recode before any model runs.
```

### Verify recoded counts against the counts a release published

[`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md):
Counts each label of `x` and compares the counts with those the source
published (a table’s margins, an annual report). When they differ and
some permutation of the labels would make them agree (up to seven labels
are searched), the labels are attached to the wrong groups, and the
message says which.

``` r

x <- rep(c("White", "Black", "Indigenous"), c(50, 30, 20))
vrfmrg(x, c(White = 50, Black = 30, Indigenous = 20))
vrfmrg(x, c(White = 30, Black = 50, Indigenous = 20), strict = FALSE)
#> vrfmrg: the counts differ 
#>       label observed published difference
#>       White       50        30         20
#>       Black       30        50        -20
#>  Indigenous       20        20          0
#> vrfmrg: recoded counts do not match the published counts. White: observed 50,
#> published 30; Black: observed 30, published 50; Indigenous: observed 20,
#> published 20. The observed counts match the published ones if the labels are
#> permuted (White -> Black, Black -> White): the labels are attached to the
#> wrong groups.
```

### Verify a categorical variable that crossed from one program to another

[`xfrvrf()`](https://rootcoder007.github.io/crimlab/reference/xfrvrf.md):
Checks a categorical column after it was exported by one program and
imported by another: the value labels that arrived must agree with the
source’s code book, and the decoded values must reproduce the counts the
source reported. When the column carries a `labels` attribute it is
decoded by it
([`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md));
otherwise by `codebook`; otherwise the values are taken as labels.

``` r

imp <- structure(c(1, 2, 2, 1, 1), labels = c(No = 1, Yes = 2))
xfrvrf(imp, c(No = 3, Yes = 2), codebook = c("1" = "No", "2" = "Yes"))
#> xfrvrf: as the source reported 
#> code book: agrees 
#> vrfmrg: the counts match 
#>  label observed published difference
#>     No        3         3          0
#>    Yes        2         2          0
xfrvrf(c("No", "Yes", "Yes"), c(No = 1, Yes = 2))
#> xfrvrf: as the source reported 
#> code book: nothing to compare 
#> vrfmrg: the counts match 
#>  label observed published difference
#>     No        1         1          0
#>    Yes        2         2          0
```

## Text

### Transliterate text to plain ASCII

[`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md):
Converts text to plain 7-bit ASCII: the common Latin accents and
ligatures are folded to their nearest ASCII letters (an accented capital
A becomes `"A"`, sharp s `"ss"`, the ae ligature `"ae"`, thorn `"th"`),
and every character with no folding is dropped. Use it for destinations
that accept only ASCII (a file name, a fixed-width export, a legacy
system), so that a name written with an accented letter arrives with the
plain letter rather than as an encoding error or a question mark.

``` r

toascii("Prof. \u00c1ngela Zorro Medina")
#> [1] "Prof. Angela Zorro Medina"
toascii(c("Se\u00e1n", "Zo\u00eb", "na\u00efve", "\u0141\u00f3d\u017a"))
#> [1] "Sean"  "Zoe"   "naive" "Lodz"
all(charToRaw(toascii("caf\u00e9")) < 128)
#> [1] TRUE
toascii(c(NA, "two\nlines"))
#> [1] ""         "twolines"
```

### Use text as it is, falling back to ASCII when it is not valid UTF-8

[`asciifb()`](https://rootcoder007.github.io/crimlab/reference/asciifb.md):
Returns `x` unchanged when it is valid UTF-8 (so an accented name is
kept), and transliterates to plain ASCII with
[`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md)
only the elements that are not valid UTF-8 (an encoding error), or every
element when `force = TRUE` (for an ASCII-only destination). Names keep
their accents wherever UTF-8 is supported and degrade gracefully where
it is not.

``` r

asciifb("\u00c1ngela")
#> [1] "Ángela"
asciifb("\u00c1ngela", force = TRUE)
#> [1] "Angela"
bad <- rawToChar(as.raw(c(0x63, 0x61, 0x66, 0xe9)))
validUTF8(bad)
#> [1] FALSE
asciifb(bad)
#> [1] "caf"
```

### Write text as UTF-8, falling back to ASCII where it is not UTF-8

[`wrtxtfb()`](https://rootcoder007.github.io/crimlab/reference/wrtxtfb.md):
Writes `text` to `path` as UTF-8, one element per line, for generated
text files (a README, a citation, a script) whose content may hold names
with accents. An element that is not valid UTF-8 is written as its ASCII
transliteration
([`asciifb()`](https://rootcoder007.github.io/crimlab/reference/asciifb.md)),
so the file is always valid UTF-8.

``` r

p <- wrtxtfb(c("line one", "line two"), tempfile(fileext = ".txt"))
readLines(p)
#> [1] "line one" "line two"
q <- wrtxtfb("Author: \u00c1ngela Zorro Medina", tempfile())
readBin(q, "raw", 20)
#>  [1] 41 75 74 68 6f 72 3a 20 c3 81 6e 67 65 6c 61 20 5a 6f 72 72
try(wrtxtfb("x", c("a.txt", "b.txt")))
#> Error : wrtxtfb: `path` must be one non-empty file path.
```
