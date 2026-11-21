# Open data from the publisher's portal

Police services, courts and corrections ministries publish through a few
kinds of portal. These functions ask a portal where a dataset's files
are now, and download one. `ckanres()` lists the files of a CKAN dataset
and `ckanls()` searches a CKAN portal; `socres()` describes a Socrata
dataset and gives the address of its CSV export; `arcres()` describes an
ArcGIS feature layer and gives the address of a query for its rows;
`srcget()` downloads an address to a file and returns its size and
SHA-256, to be recorded with the analysis; `scres()` gives the address
of a Statistics Canada table. `srcls()` lists sources the package knows
by a short key (see Details), and `srcurl()` resolves one of them.

## Usage

``` r
ckanres(portal, id)

ckanls(portal, q, rows = 50L)

socres(domain, id)

arcres(layer)

srcget(url, dest, over = FALSE, wait = 600)

srcls()

srcurl(key)

scres(pid, lang = c("en", "fr"))
```

## Arguments

- portal:

  The address of a CKAN portal, such as `"https://data.ontario.ca"`.

- id:

  A dataset's identifier on its portal: a CKAN package name or id, a
  Socrata four-by-four (`"ijzp-q8t2"`).

- q:

  Search words.

- rows:

  The most datasets `ckanls()` returns, 1 to 1000.

- domain:

  The host of a Socrata portal, such as `"data.cityofchicago.org"`.

- layer:

  The address of an ArcGIS feature layer, ending in its number
  (`.../FeatureServer/0`).

- url:

  The address to download.

- dest:

  The file to write.

- over:

  `TRUE` replaces `dest` when it exists.

- wait:

  The seconds to wait for the file (`getOption("timeout")` when that is
  longer).

- key:

  A source's key, as `srcls()` lists it.

- pid:

  The number of a Statistics Canada table: eight digits, with or without
  the hyphens of `"35-10-0177-01"` (the last two digits, the view, are
  not part of it).

- lang:

  `"en"` (the default) or `"fr"`: the language of the file.

## Value

`ckanres()` a data frame with a row per file: `name`, `format`, `url`,
`modified` and `bytes` (`NA` where the portal does not say). `ckanls()`
a data frame with a row per dataset: `id`, `title`, `files` (their
number) and `licence`. `socres()` and `arcres()` a list: `name`, `url`
(the CSV export, or the GeoJSON query), `columns` (their names),
`updated` (a `POSIXct`, or `NA`), and for a layer `rows` and `max`.
`srcget()`, invisibly, a list: `path`, `bytes`, `sha256` and `url`.
`scres()` the address of the table's zipped CSV file. `srcls()` a data
frame: `key`, `publisher`, `what`, `kind` (`"ckan"`, `"socrata"`,
`"arcgis"`, `"statcan"`, `"url"` or `"page"`), `portal` and `id`.
`srcurl()` one address: for a CKAN dataset the first CSV file's (else
the first file's), for a page or a portal its own.

## Details

Nothing is cached and nothing about a dataset is fixed in the package: a
file the publisher has replaced is the one found. The files are the
publisher's, under the licence the publisher states on its portal
(`ckanls()` gives it for a CKAN dataset); read it before redistributing
a file.

The sources of `srcls()`. In Canada: the Toronto Police Service's public
safety data (ArcGIS layers: major crime indicators, homicides,
shootings, neighbourhood crime rates, race-based data on use of force
and on arrests and strip searches, Mental Health Act apprehensions), the
City of Toronto's police annual statistical report and Montreal's
criminal acts (CKAN), Ontario's data on people in custody (the Jahn
settlement data), on use of force in correctional institutions and on
people under correctional supervision (CKAN), Statistics Canada's tables
of police-reported crime, the crime severity index, homicide, hate
crime, adult criminal courts and correctional services, the director's
reports of Ontario's Special Investigations Unit (one page each:
[`ovrget()`](https://rootcoder007.github.io/crimlab/reference/ovrget.md))
and the Vancouver Police Department's crime data, which is behind a form
that accepts its terms. In the United States: the police departments of
New York City, Chicago and Los Angeles (Socrata), the Crime Data
Explorer of the Federal Bureau of Investigation and the National Crime
Victimization Survey (pages). Elsewhere: the latest archive of the
United Kingdom's street-level crime, outcomes and stop and search (over
a gigabyte), Eurostat's police-recorded offences, and the open datasets
of the New South Wales Bureau of Crime Statistics and Research (a page).
The last rows are portals to search with `ckanls()`: Canada, Ontario,
Toronto, Montreal, the United Kingdom and Australia. A source with no
file of its own (a page, a portal) resolves to its address.

An ArcGIS layer returns at most `max` rows to one query (the layer's own
limit, 1000 or 2000 as a rule): `arcres()` reports that limit and the
layer's row count, and a larger layer is read in pages by adding
`&resultOffset=` to the address.

Errors: a `portal`, `domain`, `layer` or `url` that is not one `http://`
or `https://` address (a `domain` may be a bare host name); an `id`, `q`
or `key` that is not one non-empty string; a `pid` that is not a table
number; a portal that cannot be reached, does not answer in time or
answers with an error; an answer that is not what that kind of portal
sends; a dataset the portal does not have; a `dest` that is not one file
path, or that exists without `over = TRUE`; a `key` `srcls()` does not
list.

## See also

[`odtab()`](https://rootcoder007.github.io/crimlab/reference/odtab.md)
for the tables the project's service mirrors and checks;
[`ovrget()`](https://rootcoder007.github.io/crimlab/reference/ovrget.md)
for one oversight report.

## Examples

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
if (FALSE) { # \dontrun{
# needs the internet
f <- ckanres("https://data.ontario.ca",
             "jahn-settlement-data-inmates-ontario")
head(f[, c("name", "format", "modified")])
s <- socres("data.cityofchicago.org", "dpt3-jri9")
s$columns
g <- srcget(srcurl("tps-homicides"), tempfile(fileext = ".geojson"))
g$sha256
scres("35-10-0177-01")
} # }
try(srcurl("nowhere"))
#> Error : srcurl: `nowhere` is not a key of srcls().
```
