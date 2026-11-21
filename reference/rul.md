# Declare rules a table must satisfy, and check them

`rul()` declares one rule: a name, a function that says whether the data
pass, the column it reads (or the whole table), how serious a failure is
and the message to give. `valrules()` runs a list of rules on a data
frame and returns the ones that failed, with the rows that broke them.
The ready-made rules
([`ruleinset()`](https://rootcoder007.github.io/crimlab/reference/ruleinset.md)
and the others) cover the usual checks, so that they are written once
and not by hand.

## Usage

``` r
rul(name, expr, col = NULL, sev = c("warning", "fatal"), msg = NULL)

valrules(data, rules)
```

## Arguments

- name:

  The rule's name: one non-empty string. It names the rule in the
  result.

- expr:

  A function of one argument that returns `TRUE`/`FALSE` for the table
  or one logical value per row.

- col:

  The column the rule reads, or `NULL` (the default) for a rule on the
  whole data frame.

- sev:

  `"warning"` (the default) or `"fatal"`: how serious a failure is.
  `valrules()` records it; the caller decides what follows.

- msg:

  The message of a failure, or `NULL` for "Rule '' failed".

- data:

  The table the rules are run on: a data frame.

- rules:

  One rule or a list of rules.

## Value

`rul()` an object of class `clrule`, a list: `rule` (`"fn"`), `name`,
`col` (when given), `sev`, `msg` and `expr`. `valrules()` an object of
class `valrules`: a list with one element per failed rule, named by the
rule, each a list of `severity`, `message` and, for a row-wise rule,
`nfail` (the number of failing rows) and `rows` (the first ten of them).
A table that passes every rule gives an empty list.

## Details

The function `expr` receives the column `col` as a vector, or the whole
data frame when `col` is `NULL`. It returns either one `TRUE` or `FALSE`
for the table, or one value per row; in the row-wise form a missing
answer counts as a failure, since a rule that cannot decide has not been
met. A rule for a column the table does not have is skipped. A rule
whose function stops with an error, or returns something that is not
logical, is reported as a failure that says so: a broken rule never
passes silently. Two rules of one name leave one entry, the later one.

A ready-made rule holds no function: it is a list of its kind (`rule`)
and what it was given, and `valrules()` also takes it written by hand,
as `list(rule = "btw", col = "year", lo = 2000, hi = 2026)`. That form
is the one the other languages of the package use, where `valrules()`
and the ready-made rules exist too; a rule made by `rul()` holds an R
function and runs in R only.

Errors: a `name` or `col` that is not one non-empty string; `expr` not a
function; `sev` not `"warning"` or `"fatal"`; `data` not a data frame;
`rules` not a rule or a list of rules.

## See also

[`ruleinset()`](https://rootcoder007.github.io/crimlab/reference/ruleinset.md)
for the ready-made rules;
[`catchk()`](https://rootcoder007.github.io/crimlab/reference/catchk.md)
and
[`prechk()`](https://rootcoder007.github.io/crimlab/reference/prechk.md)
for other checks before an analysis.

## Examples

``` r
d <- data.frame(id = c(1, 2, 2, 4), age = c(34, -1, 51, NA),
                charge = c("assault", "theft", "theft", "arson"))
r <- list(
  rul("age_not_negative", function(v) is.na(v) | v >= 0, col = "age",
      msg = "an age is negative"),
  ruleuniq("id", sev = "fatal"),
  rulenull("age"),
  ruleinset("charge", c("assault", "theft", "drugs"))
)
v <- valrules(d, r)
v
#> valrules: 4 rules, 4 failed (1 fatal)
#>   warning  age_not_negative: an age is negative (1 of 4 rows); rows 2
#>   fatal  id_unique: Column 'id' has duplicate values (1 of 4 rows); rows 3
#>   warning  age_not_null: Column 'age' has missing values (1 of 4 rows); rows 4
#>   warning  charge_in_set: Column 'charge' has values outside the allowed set (1 of 4 rows); rows 4
v[["id_unique"]][["rows"]]
#> [1] 3
# a table that passes
length(valrules(d[1, ], r))
#> [1] 0
# a rule written as data
valrules(d, list(rule = "btw", col = "age", lo = 0, hi = 120))
#> valrules: 1 rules, 1 failed (0 fatal)
#>   warning  age_between: Column 'age' has values outside [0, 120] (1 of 4 rows); rows 2
try(rul("", function(v) TRUE))
#> Error : rul: `name` must be one non-empty string.
```
