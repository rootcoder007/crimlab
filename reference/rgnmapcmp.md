# Compare a recomputed region map against a published one

Matches two region maps on the unit and compares the named columns cell
by cell, reporting for each column how many cells disagree and the first
disagreement. Use it to confirm that a region map rebuilt from a
boundary file reproduces the one that was published.

## Usage

``` r
rgnmapcmp(pub, obs, unit, cols = NULL)
```

## Arguments

- pub, obs:

  The published and the recomputed map (data frames).

- unit:

  The name of the unit column, in both.

- cols:

  The columns to compare; `NULL` compares every column the two share
  apart from `unit`.

## Value

A data frame of class `rgnmapcmp`: a first row `"rows"` for units on one
side only, then one row per column, with `column`, `cells`, `mismatched`
and `first` (the first disagreement written out, `""` when none). Zero
`mismatched` throughout is agreement.

## Details

Units are matched as text with spaces trimmed; when a unit repeats, its
first row is used. Numeric columns are compared with
[`all.equal()`](https://rdrr.io/r/base/all.equal.html) at its default
tolerance (relative \\1.5 \times 10^{-8}\\), so a coordinate that
survived a round trip through text is not reported as a change;
everything else is compared as text after trimming. Two missing values
agree; a missing value against a present one does not.

What this establishes is that the recomputation reproduced the published
assignment, not that the assignment is right: the same method on the
same boundary file reproduces a definitional error perfectly. That is
what
[`rgnmapalt()`](https://rootcoder007.github.io/crimlab/reference/rgnmapalt.md)
is for.

## See also

[`rgnmapalt()`](https://rootcoder007.github.io/crimlab/reference/rgnmapalt.md),
[`rgnmapintg()`](https://rootcoder007.github.io/crimlab/reference/rgnmapintg.md).

## Examples

``` r
pub <- data.frame(inst = c("North Jail", "South Jail"),
                  cd = c("3557", "3520"))
obs <- pub
obs$cd[2] <- "3521"
rgnmapcmp(pub, obs, "inst")
#> rgnmapcmp: the maps disagree 
#>  column cells mismatched                                       first
#>    rows     2          0                                            
#>      cd     2          1 South Jail: published 3520, recomputed 3521
pub$lat <- c(44.39, 43.65)
obs2 <- data.frame(inst = c("North Jail", "Hill Jail"),
                   cd = c("3557", "3506"), lat = c(44.39 + 1e-12, 45.1))
rgnmapcmp(pub, obs2, "inst", cols = c("cd", "lat"))
#> rgnmapcmp: the maps disagree 
#>  column cells mismatched                      first
#>    rows     3          2 published only: South Jail
#>      cd     1          0                           
#>     lat     1          0                           
try(rgnmapcmp(pub, obs, "inst", cols = "lat"))
#> Error : rgnmapcmp: both frames need the compared columns; missing: lat.
```
