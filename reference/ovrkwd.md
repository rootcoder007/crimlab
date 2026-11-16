# Keyword indicators in oversight reports

Marks each report that mentions any word of each keyword group: the
signals of mental health and of race that analyses of oversight reports
tally, or any groups of one's own. A mention is a signal to read the
report, not a finding: a report may name a condition only to rule it
out.

## Usage

``` r
ovrkwd(x, vocab = "mhr", html = NA)
```

## Arguments

- x:

  The reports: pages (HTML) or texts, a character vector.

- vocab:

  The keyword groups: the name of a preset (`"mhr"`), a named list of
  character vectors (one group each), or a character vector (each word a
  group of its own).

- html:

  `NA` (the default) to read each element as HTML when it holds markup
  (a `<` followed by a letter, `/` or `!`), `TRUE` or `FALSE` to say.

## Value

A data frame with one row per report: a logical column per group (`TRUE`
when the report holds a word of it; `NA` for a report that is `NA`) and
`terms`, the words found, in the order of `vocab`, separated by `"; "`.

## Details

Words are matched whole and in any case (ASCII letters; others as
written); a word ending in `*` matches any word it starts
(`"schizophren*"` finds "schizophrenia"). The preset `"mhr"` has the
groups `mental_health` (mental health or illness, psychiatric,
schizophrenia, bipolar, depression, psychosis, suicidal, self-harm,
delusions, crisis intervention, "EDP", "emotionally disturbed", "MCIT")
and `race` (Black, African, Indigenous, First Nations, Metis, Inuit,
racialized, racial, racism, anti-Black, ethnic, South Asian).

## Examples

``` r
x <- c("The man had a history of schizophrenia; a crisis team came.",
       "An Indigenous woman was arrested.",
       "Nothing of either kind.")
k <- ovrkwd(x)
k
#>   mental_health  race        terms
#> 1          TRUE FALSE schizophren*
#> 2         FALSE  TRUE   Indigenous
#> 3         FALSE FALSE             
colSums(k[c("mental_health", "race")])
#> mental_health          race 
#>             1             1 
ovrkwd(x, list(force = c("taser", "conducted energy weapon", "baton")))
#>   force terms
#> 1 FALSE      
#> 2 FALSE      
#> 3 FALSE      
```
