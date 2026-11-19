# Write a change table as delimited text, JSON or Markdown

`yoycsv()` and `yoytsv()` write every column of the table, one row per
line, numbers with every digit (the shortest text that reads back as the
same number), text as it is (a code such as `046` stays `046`), `TRUE`
and `FALSE` as such and missing values as `na`; `meta` lines come first
as comments. In CSV a field holding a comma, a quote or a line break is
quoted, and an empty text field is written as `""` so that it stays
apart from a missing value written as empty. A TSV field cannot hold a
tab or a line break, so those become spaces. `yoyjs()` writes a JSON
object of the `meta` lines, the summary
([`yoysmry()`](https://rootcoder007.github.io/crimlab/reference/yoysmry.md))
and the rows, every digit kept. `yoymd()` writes a Markdown table for
reading, with the values rounded for display and Markdown's own
characters escaped.

## Usage

``` r
yoycsv(x, file = NULL, digits = NULL, na = "", meta = NULL, comment = "#")

yoytsv(x, file = NULL, digits = NULL, na = "", meta = NULL, comment = "#")

yoyjs(x, file = NULL, digits = NULL, pretty = TRUE, meta = NULL)

yoymd(x, file = NULL, digits = 1L, align = TRUE, meta = NULL)
```

## Arguments

- x:

  A change table
  ([`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md));
  `yoycsv()`, `yoytsv()` and `yoyjs()` take any data frame.

- file:

  A file path, or `NULL` to return the text. A name with no suffix is
  given the format's own (`.csv`, `.tsv`, `.json`, `.md`).

- digits:

  `NULL`, or the decimals of the percentage columns (0 to 15); in
  `yoymd()` the decimals of the change shown.

- na:

  The text of a missing value.

- meta:

  `NULL`, or lines that describe the table (its source, its settings),
  written as comments or as the table's caption.

- comment:

  The comment mark of the `meta` line.

- pretty:

  `TRUE` indents the JSON.

- align:

  `TRUE` pads the Markdown columns to a common width.

## Value

The text when `file` is `NULL`, otherwise the path written, invisibly.

## Details

The delimited and JSON writers change no value: `digits`, when given,
rounds the percentage columns (`pctchg`, `ppchg`, `pctlo`, `pctup`) and
nothing else, and only when asked.

## See also

[`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md),
[`yoyhtml()`](https://rootcoder007.github.io/crimlab/reference/yoyhtml.md).

## Examples

``` r
# simulated police-reported counts by offence and year, any agency
d <- data.frame(year = rep(2021:2023, 2),
                offence = rep(c("Assault", "Robbery"), each = 3),
                n = c(812, 845, 790, 120, 131, 98))
y <- yoypx(d$n, d$year, d$offence)
cat(yoycsv(y, meta = "simulated counts"))
#> # simulated counts
#> group,prd,value,prev,change,pctchg,flag,pctlo,pctup,p,padj,signif
#> Assault,2021,812,,,,no earlier period,,,,,
#> Assault,2022,845,812,33,4.064039408866995,,-5.603768879212168,14.72749053634146,0.43180618441086355,0.5279966429637909,FALSE
#> Assault,2023,790,845,-55,-6.508875739644971,,-15.259175908374855,3.136515563139186,0.18170438539856754,0.3634087707971351,FALSE
#> Robbery,2021,120,,,,no earlier period,,,,,
#> Robbery,2022,131,120,11,9.166666666666666,,-15.436293185543326,41.033274698316255,0.5279966429637909,0.5279966429637909,FALSE
#> Robbery,2023,98,131,-33,-25.19083969465649,,-43.01727965764448,-2.0593979283887176,0.03423381228214208,0.13693524912856833,FALSE
cat(yoymd(y))
#> | group   | prd  | value | previous | change | change % |     interval | note              |
#> | ------- | ---- | ----: | -------: | -----: | -------: | -----------: | ----------------- |
#> | Assault | 2021 |   812 |        — |      — |        — |              | no earlier period |
#> | Assault | 2022 |   845 |      812 |     33 |    +4.1% |  [-6%, +15%] |                   |
#> | Assault | 2023 |   790 |      845 |    -55 |    -6.5% |  [-15%, +3%] |                   |
#> | Robbery | 2021 |   120 |        — |      — |        — |              | no earlier period |
#> | Robbery | 2022 |   131 |      120 |     11 |    +9.2% | [-15%, +41%] |                   |
#> | Robbery | 2023 |    98 |      131 |    -33 |   -25.2% |  [-43%, -2%] |                   |
# simulated admissions to custody, with a code that must stay "046"
a <- data.frame(code = c("046", "112"), n = c(1 / 3, 2))
cat(yoytsv(a))
#> code n
#> 046  0.3333333333333333
#> 112  2
cat(yoyjs(y, pretty = FALSE))
#> {"summary":[{"group":"Assault","from":2021,"to":2023,"first":812,"last":790,"totalpct":-2.7093596059113323,"cagrpct":-1.3639820379549805,"periods":3,"up":1,"down":1,"flat":0},{"group":"Robbery","from":2021,"to":2023,"first":120,"last":98,"totalpct":-18.333333333333336,"cagrpct":-9.63038858849361,"periods":3,"up":1,"down":1,"flat":0}],"rows":[{"group":"Assault","prd":2021,"value":812,"prev":null,"change":null,"pctchg":null,"flag":"no earlier period","pctlo":null,"pctup":null,"p":null,"padj":null,"signif":null},{"group":"Assault","prd":2022,"value":845,"prev":812,"change":33,"pctchg":4.064039408866995,"flag":null,"pctlo":-5.603768879212168,"pctup":14.72749053634146,"p":0.43180618441086355,"padj":0.5279966429637909,"signif":false},{"group":"Assault","prd":2023,"value":790,"prev":845,"change":-55,"pctchg":-6.508875739644971,"flag":null,"pctlo":-15.259175908374855,"pctup":3.136515563139186,"p":0.18170438539856754,"padj":0.3634087707971351,"signif":false},{"group":"Robbery","prd":2021,"value":120,"prev":null,"change":null,"pctchg":null,"flag":"no earlier period","pctlo":null,"pctup":null,"p":null,"padj":null,"signif":null},{"group":"Robbery","prd":2022,"value":131,"prev":120,"change":11,"pctchg":9.166666666666666,"flag":null,"pctlo":-15.436293185543326,"pctup":41.033274698316255,"p":0.5279966429637909,"padj":0.5279966429637909,"signif":false},{"group":"Robbery","prd":2023,"value":98,"prev":131,"change":-33,"pctchg":-25.19083969465649,"flag":null,"pctlo":-43.01727965764448,"pctup":-2.0593979283887176,"p":0.03423381228214208,"padj":0.13693524912856833,"signif":false}]}
```
