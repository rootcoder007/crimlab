# Police oversight reports and open data

Reading the text of police oversight reports into fields, and loading
tables of the project’s open data service.

Each section is one help page: what the function is for, the example of
its help page as it runs, and the plot of its result where the package
draws one. [`?name`](https://rdrr.io/r/base/name.html) gives the
arguments and the value in full.

## Police oversight reports

### Text of a police oversight report

[`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md):
Turns the page of a report published by a police oversight body (a
civilian investigation unit, an ombudsman, an inspector) into the plain
text its other tools read: markup removed, each paragraph, heading, list
item and table row on a line of its own, so that section headings stand
alone, and character entities decoded. Use it to read or search saved
reports without a browser, or before the rules that code them.

``` r

ovrtxt("<h2>The Team</h2><p>Number of investigators assigned: 3</p>")
#> [1] "The Team\n\nNumber of investigators assigned: 3"
ovrtxt("<p>WO\n#1 Interviewed</p><p>WO&nbsp;#2 Interviewed</p>")
#> [1] "WO #1 Interviewed\n\nWO #2 Interviewed"
ovrtxt("a  b\r\nc", html = FALSE)
#> [1] "a b\nc"
```

### Dates written out in reports, as ISO dates

[`ovrdte()`](https://rootcoder007.github.io/crimlab/reference/ovrdte.md):
Finds the first date written in each text, in English or French, and
gives it as `YYYY-MM-DD`, the form that sorts and that
[`as.Date()`](https://rdrr.io/r/base/as.Date.html) reads.
`"January 5, 2023"`, `"Jan. 5th 2023"`, `"5 January 2023"`,
`"5th of January, 2023"`, `"5 janvier 2023"` and `"2023-01-05"` all give
`"2023-01-05"`.

``` r

ovrdte(c("January 5, 2023", "le 1er février 2021", "no date"))
#> [1] "2023-01-05" "2021-02-01" NA
ovrdte("On February 30, 2023 or rather March 2, 2023, the SIU ...")
#> [1] "2023-03-02"
as.Date(ovrdte("Date: April 28, 2023"))
#> [1] "2023-04-28"
```

### Field schemas for coding police oversight reports

[`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md):
The fields a preset reads from each report, and the rule, cues and
section of each: a data frame that
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
takes as `schema`. Copy it and change rows, or write one of the same
columns, to code the reports of any oversight body; nothing in the rules
is tied to one agency.

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
rbind(s, data.frame(name = "date_received", rule = "date",
                    cue = "received the complaint", alt = "", sect = "",
                    desc = "the date the complaint came in"))
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
#> 18                 date_received   date
#>                                                                                           cue
#> 1                                                                                      Case #
#> 2                                                     Police Service|Provincial Police|Police
#> 3                                                                                            
#> 4                       notified the SIU|contacted the SIU|SIU was notified|SIU was contacted
#> 5                                                                                       Date:
#> 6                                                        Number of SIU Investigators assigned
#> 7                                               Number of SIU Forensic Investigators assigned
#> 8                                                                                          WO
#> 9                                                                                          CW
#> 10                                                                                         SO
#> 11                                                                          year-old|year old
#> 12                                                                          year-old|year old
#> 13                             charged with|charges have been laid|criminal charges were laid
#> 14                                                  Director|Acting Director|Interim Director
#> 15   in the City of|in the Town of|in the Township of|in the Municipality of|in the Region of
#> 16 fracture|fractures|fractured|broken|laceration|lacerations|gunshot|stab|stabbed|burn|burns
#> 17                                                                                    Section
#> 18                                                                     received the complaint
#>                                                                       alt
#> 1                                                                        
#> 2                                                                        
#> 3                                                                        
#> 4                                                                        
#> 5                                                                        
#> 6                                                                        
#> 7                                                                        
#> 8                                                                        
#> 9                                                                        
#> 10                                                                       
#> 11                                                                       
#> 12                                                                       
#> 13 no reasonable grounds|No charges|no charges|not be laid|decline to lay
#> 14                                                                       
#> 15                                                                       
#> 16                                                                       
#> 17                                                                       
#> 18                                                                       
#>                                                                 sect
#> 1                                                                   
#> 2                                                                   
#> 3                                                 Incident Narrative
#> 4                                                                   
#> 5                                                                   
#> 6                                                                   
#> 7                                                                   
#> 8                                                                   
#> 9                                                                   
#> 10                                                                  
#> 11                                                                  
#> 12                                                                  
#> 13 Analysis and Director's Decision|Analysis and Director’s Decision
#> 14                                                                  
#> 15                                                                  
#> 16       Nature of Injuries / Treatment|Nature of Injuries/Treatment
#> 17                                              Relevant Legislation
#> 18                                                                  
#>                                                    desc
#> 1                         the case number of the report
#> 2  the police service whose officials were investigated
#> 3                              the date of the incident
#> 4                        the date the unit was notified
#> 5                   the date of the director's decision
#> 6                  the number of investigators assigned
#> 7         the number of forensic investigators assigned
#> 8          the number of witness officials (WO #1, ...)
#> 9         the number of civilian witnesses (CW #1, ...)
#> 10         the number of subject officials (SO #1, ...)
#> 11                       the age of the affected person
#> 12             the sex or gender of the affected person
#> 13       whether the director found grounds for charges
#> 14                 the name of the director who decided
#> 15         the municipality where the incident happened
#> 16                 the injuries the affected person had
#> 17               the provisions the director considered
#> 18                       the date the complaint came in
```

### Code police oversight reports into fields

[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md):
Reads the fields of a schema from each report: the dates, the counts of
officials and witnesses, the service, the outcome and whatever else the
schema names, each by its rule (see
[`ovrsch()`](https://rootcoder007.github.io/crimlab/reference/ovrsch.md)).
The rules are deterministic and offline, so the same report always gives
the same row, and a folder of saved reports becomes a data set whose
every value can be traced to the words that gave it.

``` r

r <- paste(
  "<h1>Director's Report - Case # 23-OFD-001</h1>",
  "<p>On January 6, 2023, the Barrie Police Service contacted the SIU.",
  "In the night, a 34-year-old man was injured.</p>",
  "<p>Number of SIU Investigators assigned: 3</p>",
  "<p>WO #1 Interviewed</p><p>WO #2 Interviewed</p>",
  "<p>Joseph Martino</p><p>Director</p>"
)
p <- ovrpar(r)
p[, c("case_number", "police_service", "date_siu_notified_iso",
      "siu_investigators", "number_of_witness_officials",
      "age_affected", "directors_name")]
#>   case_number        police_service date_siu_notified_iso siu_investigators
#> 1  23-OFD-001 Barrie Police Service            2023-01-06                 3
#>   number_of_witness_officials age_affected directors_name
#> 1                           2           34 Joseph Martino
```

### The number of subject officials in an oversight report

[`ovrsoc()`](https://rootcoder007.github.io/crimlab/reference/ovrsoc.md):
Counts the subject officials of each report, the officers whose conduct
the oversight body investigated, and says which rule gave the count and
on what words, so that the number can be checked against the report. Use
it where a report has no roster
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
can read, or to check a count coded by hand.

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

### Keyword indicators in oversight reports

[`ovrkwd()`](https://rootcoder007.github.io/crimlab/reference/ovrkwd.md):
Marks each report that mentions any word of each keyword group: the
signals of mental health and of race that analyses of oversight reports
tally, or any groups of one’s own. A mention is a signal to read the
report, not a finding: a report may name a condition only to rule it
out.

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

### Check a coded table of oversight reports

[`ovrchk()`](https://rootcoder007.github.io/crimlab/reference/ovrchk.md):
Flags the cells of a table coded from oversight reports (by
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md),
by hand, or by any other means) that are not of the form their field
calls for, so that a parsing slip or a typing error is found before the
table is analysed. Nothing is changed.

``` r

d <- data.frame(case_number = c("23-OFD-001", "23-OFD-1", NA),
                date_of_incident_iso = c("2023-01-05", "2023-02-30", ""),
                number_of_subject_officials = c(1, 2.5, NA),
                charges_recommended = c(FALSE, NA, TRUE))
ovrchk(d, shape = c(case_number = "99-AAA-999|99-AA-999"))
#>   row issues
#> 1   1      0
#> 2   2      3
#> 3   3      0
#>                                                                                 what
#> 1                                                                                   
#> 2 date_of_incident_iso: date; number_of_subject_officials: count; case_number: shape
#> 3
```

## Open data

### Tables of the open data service

[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md),
[`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md):
[`odls()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
lists the tables of the project’s open data service and
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
loads one by its key, as the service serves it at the time of the call.
Everything about a table (its address, rows, columns and SHA-256) is
read from the service’s manifest at run time, so a table that its
publisher has replaced is loaded as it is now, with no new release of
crimlab.

``` r

## Not run:
##D # needs the internet
##D tabs <- odls()
##D head(tabs[, c("key", "rows")])
##D x <- odtab(tabs$key[!is.na(tabs$sha256)][1L])
##D str(x)
## End(Not run)
```
