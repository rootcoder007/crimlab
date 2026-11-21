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
s2 <- rbind(s, data.frame(name = "date_received", rule = "date",
                          cue = "received the complaint", alt = "",
                          sect = "",
                          desc = "the date the complaint came in"))
tail(s2[, c("name", "rule", "cue")], 2)
#>                    name  rule                    cue
#> 17 relevant_legislation lines                Section
#> 18        date_received  date received the complaint
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

### The tables of the data service about Indigenous peoples

[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md):
Lists the tables of the project’s open data service drawn from
Statistics Canada’s tables on Indigenous peoples and from the 2021
Census Indigenous Population Profile, with their descriptions, sizes and
licences, optionally those whose description matches a pattern.
[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
loads one by its key.

``` r

## Not run:
##D # needs the internet
##D t <- indtabs("mental health")
##D t[, c("key", "rows")]
## End(Not run)
try(indtabs(pat = ""))
#> Error : indtabs: `pat` must be one non-empty string.
```

### Statistics Canada tables on Indigenous peoples, by subject

[`indiex()`](https://rootcoder007.github.io/crimlab/reference/indiex.md):
The tables of the project’s open data service that hold what is most
often asked for about First Nations people, Metis and Inuit: health
indicators, perceived mental health, wellbeing, suicidal thoughts, the
census counts of population and of households, life expectancy and
infant mortality.
[`indiex()`](https://rootcoder007.github.io/crimlab/reference/indiex.md)
finds the tables of one subject among those
[`indtabs()`](https://rootcoder007.github.io/crimlab/reference/indtabs.md)
lists, and with `load = TRUE` loads the first.

``` r

## Not run:
##D # needs the internet
##D indiex("mentalhealth")[, c("key", "rows")]
##D x <- indiex("lifeexp", load = TRUE)
## End(Not run)
try(indiex("weather"))
#> Error : indiex: `name` must be one of "health", "mentalhealth", "wellbeing", "suicide", "censuspop", "censushh", "lifeexp", "infantmort", not "weather".
```

### Open data from the publisher’s portal

[`ckanres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md),
[`ckanls()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md),
[`socres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md),
[`arcres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
and others: Police services, courts and corrections ministries publish
through a few kinds of portal. These functions ask a portal where a
dataset’s files are now, and download one.
[`ckanres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
lists the files of a CKAN dataset and
[`ckanls()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
searches a CKAN portal;
[`socres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
describes a Socrata dataset and gives the address of its CSV export;
[`arcres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
describes an ArcGIS feature layer and gives the address of a query for
its rows;
[`srcget()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
downloads an address to a file and returns its size and SHA-256, to be
recorded with the analysis;
[`scres()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
gives the address of a Statistics Canada table.
[`srcls()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
lists sources the package knows by a short key (see Details), and
[`srcurl()`](https://rootcoder007.github.io/crimlab/reference/ckanres.md)
resolves one of them.

``` r

srcls()[, c("key", "publisher", "kind")]
#>                    key                                   publisher    kind
#> 1              tps-mci                      Toronto Police Service  arcgis
#> 2        tps-homicides                      Toronto Police Service  arcgis
#> 3        tps-shootings                      Toronto Police Service  arcgis
#> 4            tps-rates                      Toronto Police Service  arcgis
#> 5              tps-uof                      Toronto Police Service  arcgis
#> 6          tps-arrests                      Toronto Police Service  arcgis
#> 7              tps-mha                      Toronto Police Service  arcgis
#> 8   toronto-asr-crimes                             City of Toronto    ckan
#> 9  toronto-asr-arrests                             City of Toronto    ckan
#> 10     montreal-crimes                           Ville de Montreal    ckan
#> 11         siu-reports      Special Investigations Unit of Ontario    page
#> 12           vpd-crime                 Vancouver Police Department    page
#> 13        ontario-jahn   Ontario Ministry of the Solicitor General    ckan
#> 14         ontario-uof   Ontario Ministry of the Solicitor General    ckan
#> 15 ontario-supervision   Ontario Ministry of the Solicitor General    ckan
#> 16       statcan-crime                           Statistics Canada statcan
#> 17         statcan-csi                           Statistics Canada statcan
#> 18    statcan-homicide                           Statistics Canada statcan
#> 19        statcan-hate                           Statistics Canada statcan
#> 20      statcan-courts                           Statistics Canada statcan
#> 21  statcan-admissions                           Statistics Canada statcan
#> 22     statcan-custody                           Statistics Canada statcan
#> 23     nypd-complaints             New York City Police Department socrata
#> 24 nypd-complaints-ytd             New York City Police Department socrata
#> 25        nypd-arrests             New York City Police Department socrata
#> 26    nypd-arrests-ytd             New York City Police Department socrata
#> 27           nypd-hate             New York City Police Department socrata
#> 28          cpd-crimes                   Chicago Police Department socrata
#> 29         cpd-arrests                   Chicago Police Department socrata
#> 30         lapd-crimes               Los Angeles Police Department socrata
#> 31             fbi-cde             Federal Bureau of Investigation    page
#> 32            bjs-ncvs                Bureau of Justice Statistics    page
#> 33           uk-police                Home Office (data.police.uk)     url
#> 34      eurostat-crime                                    Eurostat     url
#> 35          nsw-bocsar NSW Bureau of Crime Statistics and Research    page
#> 36       portal-canada                        Government of Canada    ckan
#> 37      portal-ontario                       Government of Ontario    ckan
#> 38      portal-toronto                             City of Toronto    ckan
#> 39     portal-montreal                           Ville de Montreal    ckan
#> 40           portal-uk            Government of the United Kingdom    ckan
#> 41    portal-australia                       Australian Government    ckan
## Not run:
##D # needs the internet
##D f <- ckanres("https://data.ontario.ca",
##D              "jahn-settlement-data-inmates-ontario")
##D head(f[, c("name", "format", "modified")])
##D s <- socres("data.cityofchicago.org", "dpt3-jri9")
##D s$columns
##D g <- srcget(srcurl("tps-homicides"), tempfile(fileext = ".geojson"))
##D g$sha256
##D scres("35-10-0177-01")
## End(Not run)
try(srcurl("nowhere"))
#> Error : srcurl: `nowhere` is not a key of srcls().
```

### Fetch one police oversight report

[`ovrget()`](https://rootcoder007.github.io/crimlab/reference/ovrget.md):
Downloads a director’s report of Ontario’s Special Investigations Unit
by its report number, as the page the Unit serves, to a file.
[`ovrtxt()`](https://rootcoder007.github.io/crimlab/reference/ovrtxt.md)
then gives its text and
[`ovrpar()`](https://rootcoder007.github.io/crimlab/reference/ovrpar.md)
codes it into fields.

``` r

## Not run:
##D # needs the internet
##D f <- ovrget(1234, tempfile(fileext = ".html"))
##D substr(ovrtxt(readLines(f$path, warn = FALSE)), 1, 200)
## End(Not run)
try(ovrget(0, tempfile()))
#> Error : ovrget: `drid` must be one positive report number (the drid of the report's address).
```
