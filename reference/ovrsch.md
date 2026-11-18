# Field schemas for coding police oversight reports

The fields a preset reads from each report, and the rule, cues and
section of each: a data frame that
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
takes as `schema`. Copy it and change rows, or write one of the same
columns, to code the reports of any oversight body; nothing in the rules
is tied to one agency.

## Usage

``` r
ovrsch(preset = "siu")
```

## Arguments

- preset:

  The name of a preset: `"siu"`.

## Value

A data frame with one row per field: `name`, `rule`, `cue`, `alt`,
`sect` and `desc` (what the field records).

## Details

A field is read by one rule. `cue` holds the phrases the rule looks for,
alternatives separated by `"|"`, matched as written (case counts) and as
whole words at their ends. With `sect` (headings separated by `"|"`),
the rule reads only the section under the last line that is one of those
headings (a page that opens with a table of contents repeats each
heading), up to the next line that looks like a heading (a short line
starting with a capital, with no digit, not ending in a stop or a
comma); a report without the heading leaves the field missing. The
rules:

- `after`: the text after a cue (past spaces and a colon), up to a
  comma, a semicolon, a full stop or the end of the line.

- `num`: the first number (digits, or a number word in English or French
  up to twenty) within 40 characters after a cue on its line.

- `date`: the first date
  ([`ovrdte()`](https://rootcoder007.github.io/crimlab/reference/ovrdte.md))
  in the sentence of a cue, or with no cue the first date in the
  section.

- `tag`: the highest number of a numbered tag (`"WO #3"`, `"WO 3"`,
  `"AT no 3"`): the count of the people the tags number.

- `sent`: the sentence that holds the first cue.

- `prev`: the line before the last line that is a cue (a signature: the
  name above a title).

- `lines`: every line that starts with a cue, distinct, joined by
  `"; "`.

- `name`: a name of up to four capitalised words that ends in a cue
  (`"Barrie Police Service"` for the cue `"Police Service"`), without a
  leading "The".

- `flag`: `FALSE` when a phrase of `alt` occurs, else `TRUE` when a cue
  does.

- `before`: the number just before a cue (`"34"` in `"34-year-old"`).

- `word`: the word just after a cue (`"man"` in `"34-year-old man"`).

The preset `"siu"` reads the director's reports of Ontario's Special
Investigations Unit as published in English since 2020: the case number,
the police service, the dates of the incident, of the notification and
of the decision, the investigators assigned, the witness officials,
civilian witnesses and subject officials, the age and sex of the
affected person, whether the director found grounds for charges, the
director's name, the municipality, the injuries and the provisions
considered.

## Examples

``` r
s <- ovrsch()
s[, c("name", "rule")]
#>                             name   rule
#> 1                    case_number  after
#> 2                 police_service   name
#> 3           date_of_incident_iso   date
#> 4          date_siu_notified_iso   date
#> 5  date_of_director_decision_iso   date
#> 6              siu_investigators    num
#> 7    siu_forensics_investigators    num
#> 8    number_of_witness_officials    tag
#> 9   number_of_civilian_witnesses    tag
#> 10   number_of_subject_officials    tag
#> 11                  age_affected before
#> 12           sex_gender_affected   word
#> 13           charges_recommended   flag
#> 14                directors_name   prev
#> 15              location_of_call  after
#> 16             specific_injuries   sent
#> 17          relevant_legislation  lines
# a field of one's own: the date a complaint came in
s2 <- rbind(s, data.frame(name = "date_received", rule = "date",
                          cue = "received the complaint", alt = "",
                          sect = "",
                          desc = "the date the complaint came in"))
tail(s2[, c("name", "rule", "cue")], 2)
#>                    name  rule                    cue
#> 17 relevant_legislation lines                Section
#> 18        date_received  date received the complaint
```
