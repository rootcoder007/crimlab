# Text of a police oversight report

Turns the page of a report published by a police oversight body (a
civilian investigation unit, an ombudsman, an inspector) into the plain
text its other tools read: markup removed, each paragraph, heading, list
item and table row on a line of its own, so that section headings stand
alone, and character entities decoded. Use it to read or search saved
reports without a browser, or before the rules that code them.

## Usage

``` r
ovrtxt(x, html = NA)
```

## Arguments

- x:

  The pages or texts: a character vector. `NA` stays `NA`.

- html:

  `NA` (the default) to read each element as HTML when it holds markup
  (a `<` followed by a letter, `/` or `!`), `TRUE` or `FALSE` to say.

## Value

A character vector the length of `x`, in UTF-8.

## Details

Scripts, styles and comments are dropped. The start and end of a block
element (a paragraph, a heading, a division, a list item, a table row, a
line break and the like) end a line; every other tag is a space, and the
page's own line breaks are spaces, as a browser shows them. Named
entities common on report pages and every numeric entity are decoded
(curly quotes and dashes to their plain forms). Spaces are collapsed,
lines trimmed, and at most one empty line is kept between paragraphs.
Plain text (with `html = FALSE`, or with no markup when `html = NA`) is
only tidied that way, its line breaks kept. Nothing is fetched and
nothing is written; the pass is linear in the length of the page.

## Examples

``` r
ovrtxt("<h2>The Team</h2><p>Number of investigators assigned: 3</p>")
#> [1] "The Team\n\nNumber of investigators assigned: 3"
ovrtxt("<p>WO\n#1 Interviewed</p><p>WO&nbsp;#2 Interviewed</p>")
#> [1] "WO #1 Interviewed\n\nWO #2 Interviewed"
ovrtxt("a  b\r\nc", html = FALSE)
#> [1] "a b\nc"
```
