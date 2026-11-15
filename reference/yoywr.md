# Write a change table in the format its file name says

Picks the writer from the extension of `file` (`.html` or `.htm`,
`.pdf`, `.csv`, `.tsv` or `.tab`, `.json`, `.md` or `.markdown`), or
from `format`, and passes the other arguments to it.

## Usage

``` r
yoywr(x, file, format = "auto", ...)
```

## Arguments

- x:

  A change table
  ([`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)).

- file:

  The file path.

- format:

  `"auto"` (from the extension), or one of `"html"`, `"pdf"`, `"csv"`,
  `"tsv"`, `"json"`, `"markdown"`.

- ...:

  The arguments of the chosen writer.

## Value

`file`, invisibly.

## See also

[`yoycsv()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md),
[`yoyhtml()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md).

## Examples

``` r
# simulated police-reported counts by offence and year
y <- yoypx(c(812, 845, 790), 2021:2023)
f <- tempfile(fileext = ".csv")
yoywr(y, f)
readLines(f)
#> [1] "prd,value,prev,change,pctchg,flag,pctlo,pctup,p,padj,signif"                                                            
#> [2] "2021,812,,,,no earlier period,,,,,"                                                                                     
#> [3] "2022,845,812,33,4.064039408866995,,-5.603768879212168,14.72749053634146,0.43180618441086355,0.43180618441086355,FALSE"  
#> [4] "2023,790,845,-55,-6.508875739644971,,-15.259175908374855,3.136515563139186,0.18170438539856754,0.3634087707971351,FALSE"
unlink(f)
```
