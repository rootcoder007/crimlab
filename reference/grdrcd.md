# Recode categories by name

`grdrcd()` recodes a categorical column by an explicit mapping from old
to new labels and refuses any observed value the mapping does not name:
silent pass-through is how group labels get corrupted. `deccodes()`
decodes numeric category codes with a code-to-label dictionary (codes as
names), the step that turns the codes of an imported file into its
categories. `relbl()` relabels a variable by name, never by position,
and returns a factor whose levels follow the mapping.

## Usage

``` r
grdrcd(x, map, keep = character())

deccodes(x, vlabs, keepna = TRUE)

relbl(x, map, keep = character())
```

## Arguments

- x:

  The labels or codes: a vector or a factor.

- map, vlabs:

  A named character vector: old label (or code) = new label.

- keep:

  Values that pass through unchanged.

- keepna:

  `FALSE` refuses missing codes.

## Value

`grdrcd()` and `deccodes()`: a character vector the length of `x`.
`relbl()`: a factor.

## Details

Assigning labels by position is how category codes get rotated: when the
labels are listed in alphabetical order and the codes are not, every
group receives its neighbour's label. Mapping by name makes the
assignment independent of any order.

In `grdrcd()` a value listed in `keep` passes through unchanged and
missing values stay missing. In `relbl()` every level (the levels of a
factor, or the sorted distinct values) must be mapped or kept; the
levels of the result are the new labels in the order of the mapping,
then the kept values that occur. The result carries an attribute
`recode`: the mapping and the kept values.

## See also

[`grdlvls()`](https://rootcoder007.github.io/crimlab/reference/grdlvls.md),
[`declbld()`](https://rootcoder007.github.io/crimlab/reference/declbld.md),
[`vrfrcd()`](https://rootcoder007.github.io/crimlab/reference/vrfrcd.md),
[`vrfmrg()`](https://rootcoder007.github.io/crimlab/reference/vrfmrg.md).

## Examples

``` r
# simulated race codes of an imported file
codes <- c(1, 2, 2, 3, 1, NA)
deccodes(codes, c("1" = "White", "2" = "Black", "3" = "Indigenous"))
#> [1] "White"      "Black"      "Black"      "Indigenous" "White"     
#> [6] NA          
#> attr(,"recode")
#> attr(,"recode")$mapping
#>            1            2            3 
#>      "White"      "Black" "Indigenous" 
#> 
#> attr(,"recode")$kept
#> character(0)
#> 
grdrcd(c("M", "F", "F", "X"), c(M = "Male", F = "Female"), keep = "X")
#> [1] "Male"   "Female" "Female" "X"     
#> attr(,"recode")
#> attr(,"recode")$mapping
#>        M        F 
#>   "Male" "Female" 
#> 
#> attr(,"recode")$kept
#> [1] "X"
#> 
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
