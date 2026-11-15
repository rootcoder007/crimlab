# Render a change table to HTML or PDF

`yoyhtml()` writes one self-contained HTML page (every style inlined,
nothing fetched, so the file renders the same years later): the table, a
card per series with its change over the span
([`yoysmry()`](https://rootcoder007.github.io/crimlab/reference/yoysmry.md))
and the notes. `yoypdf()` writes a PDF of the table, paginated when it
does not fit, with the standard Helvetica font and no other program:
crimlab writes the file itself.

## Usage

``` r
yoyhtml(
  x,
  file,
  title = "Year-over-year change",
  sub = "",
  pal = "diverging",
  better = "neutral",
  digits = 1L,
  bars = TRUE,
  intv = TRUE,
  notes = NULL,
  lag = 1
)

yoypdf(
  x,
  file,
  title = "Year-over-year change",
  sub = "",
  pal = "diverging",
  better = "neutral",
  digits = 1L,
  intv = TRUE,
  notes = NULL,
  width = 11,
  height = 8.5,
  lag = 1
)
```

## Arguments

- x:

  A change table
  ([`yoypx()`](https://rootcoder007.github.io/crimlab/reference/yoypx.md)).

- file:

  The file path.

- title:

  The heading (and the page title).

- sub:

  The line under the heading.

- pal:

  The palette: `"diverging"`, `"safe"` or `"mono"`.

- better:

  Which direction is better: `"neutral"`, `"higher"` or `"lower"`.

- digits:

  The decimals of the change shown.

- bars:

  `TRUE` draws a bar under each change, proportional to its size.

- intv:

  `TRUE` shows the interval column of a count table.

- notes:

  `NULL`, or footnotes.

- lag:

  The periods per step of the table, for the summary and the note on the
  first periods.

- width, height:

  The PDF page in inches, each from 2 to 200.

## Value

`file`, invisibly.

## Details

Colour encodes the direction of an improvement, not the sign of the
change: with `better = "lower"` (fewer assaults, fewer days in
segregation) a rise is drawn in the colour for worse, with `"higher"` (a
completion rate) in the colour for better, and with `"neutral"` (the
default) every change in the neutral colour. The palette `"safe"` (blue
and orange) stays apart in the common forms of colour blindness;
`"mono"` is for print. The first period of each series, with nothing to
compare with, is said once under the table rather than on every first
row. The PDF holds ASCII only: dashes and arrows become `-` and `+`,
other characters their transliteration
([`toascii()`](https://rootcoder007.github.io/crimlab/reference/toascii.md)).

## See also

[`yoywr()`](https://rootcoder007.github.io/crimlab/reference/yoywr.md),
[`yoycsv()`](https://rootcoder007.github.io/crimlab/reference/yoycsv.md).

## Examples

``` r
# simulated police-reported counts by offence and year, any agency
d <- data.frame(year = rep(2021:2023, 2),
                offence = rep(c("Assault", "Robbery"), each = 3),
                n = c(812, 845, 790, 120, 131, 98))
y <- yoypx(d$n, d$year, d$offence)
f <- tempfile(fileext = ".html")
yoyhtml(y, f, title = "Reported offences (simulated)", better = "lower")
g <- tempfile(fileext = ".pdf")
yoypdf(y, g, title = "Reported offences (simulated)")
# simulated days in segregation by fiscal year
s <- yoylbl(yoypx(c(4100, 3900, 4400), 2021:2023), fylbl)
yoyhtml(s, f, title = "Days in segregation (simulated)", better = "lower")
unlink(c(f, g))
```
