# The number of subject officials in an oversight report

Counts the subject officials of each report, the officers whose conduct
the oversight body investigated, and says which rule gave the count and
on what words, so that the number can be checked against the report. Use
it where a report has no roster
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
can read, or to check a count coded by hand.

## Usage

``` r
ovrsoc(
  x,
  tag = "SO",
  term = c("subject official", "subject officer"),
  wtag = "WO",
  head = c("Subject Officials", "Subject Officers"),
  html = NA
)
```

## Arguments

- x:

  The reports: pages (HTML) or texts, a character vector.

- tag:

  The tag of a subject official in the report: one text.

- term:

  The words for a subject official, singular: text.

- wtag:

  The tag of a witness official: one text.

- head:

  The headings of the roster of subject officials: text.

- html:

  `NA` (the default) to read each element as HTML when it holds markup
  (a `<` followed by a letter, `/` or `!`), `TRUE` or `FALSE` to say.

## Value

A data frame with one row per report: `count` (whole number, `NA` when
unresolved), `rule` (1 to 6, as above) and `reason` (the rule and the
words it read).

## Details

The rules run in order and the first that applies gives the count:

1.  the roster under a heading of `head`: the highest numbered tag there
    (`"SO #2"`);

2.  the highest numbered tag in the whole report, trusted only when the
    tag numbered 1 is there too;

3.  a number before a plural of `term` (`"two subject officials"`,
    `"3 subject officers"`);

4.  a mention of `term` in the singular, not after "no" or "not a": 1;

5.  a statement that there was none (`"no subject official"`), or
    numbered witness tags (`wtag`) with no subject tag at all: 0, a real
    answer;

6.  otherwise unresolved: `NA`.

Tags and headings are matched as written; `term` in any case. The
defaults are the words of Canadian oversight bodies; give a body's own
(`tag = "AI"`, `term = "agent impliqué"` for a French report).

## Examples

``` r
ovrsoc(c("Subject Officials\nSO #1 Interviewed\nSO #2 Declined",
         "Two subject officials were designated.",
         "Witness Officials\nWO #1\nWO #2",
         "No officers are named here."))
#>   count rule                                      reason
#> 1     2    1 the roster under the heading: highest SO #2
#> 2     2    3           a plural: 'Two subject officials'
#> 3     0    5               witness officials only: WO #2
#> 4    NA    6                                  unresolved
```
