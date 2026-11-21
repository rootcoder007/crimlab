# Data straight from its publisher. odtab() loads the tables the project's
# service mirrors and checks; these ask the publisher's own portal where a
# dataset is now (CKAN, Socrata, ArcGIS) and download it, for data the
# service does not hold or a release newer than its copy. Nothing about a
# dataset is written here but where to ask: its files, their addresses and
# their dates are read from the portal at every call.

#' Open data from the publisher's portal
#'
#' Police services, courts and corrections ministries publish through a
#' few kinds of portal. These functions ask a portal where a dataset's
#' files are now, and download one. `ckanres()` lists the files of a CKAN
#' dataset and `ckanls()` searches a CKAN portal; `socres()` describes a
#' Socrata dataset and gives the address of its CSV export; `arcres()`
#' describes an ArcGIS feature layer and gives the address of a query for
#' its rows; `srcget()` downloads an address to a file and returns its
#' size and SHA-256, to be recorded with the analysis; `scres()` gives the
#' address of a Statistics Canada table. `srcls()` lists sources the
#' package knows by a short key (see Details), and `srcurl()` resolves one
#' of them.
#'
#' Nothing is cached and nothing about a dataset is fixed in the package:
#' a file the publisher has replaced is the one found. The files are the
#' publisher's, under the licence the publisher states on its portal
#' (`ckanls()` gives it for a CKAN dataset); read it before redistributing
#' a file.
#'
#' The sources of `srcls()`. In Canada: the Toronto Police Service's
#' public safety data (ArcGIS layers: major crime indicators, homicides,
#' shootings, neighbourhood crime rates, race-based data on use of force
#' and on arrests and strip searches, Mental Health Act apprehensions),
#' the City of Toronto's police annual statistical report and Montreal's
#' criminal acts (CKAN), Ontario's data on people in custody (the Jahn
#' settlement data), on use of force in correctional institutions and on
#' people under correctional supervision (CKAN), Statistics Canada's
#' tables of police-reported crime, the crime severity index, homicide,
#' hate crime, adult criminal courts and correctional services, the
#' director's reports of Ontario's Special Investigations Unit (one page
#' each: [ovrget()]) and the Vancouver Police Department's crime data,
#' which is behind a form that accepts its terms. In the United States:
#' the police departments of New York City, Chicago and Los Angeles
#' (Socrata), the Crime Data Explorer of the Federal Bureau of
#' Investigation and the National Crime Victimization Survey (pages).
#' Elsewhere: the latest archive of the United Kingdom's street-level
#' crime, outcomes and stop and search (over a gigabyte), Eurostat's
#' police-recorded offences, and the open datasets of the New South Wales
#' Bureau of Crime Statistics and Research (a page). The last rows are
#' portals to search with `ckanls()`: Canada, Ontario, Toronto, Montreal,
#' the United Kingdom and Australia. A source with no file of its own (a
#' page, a portal) resolves to its address.
#'
#' An ArcGIS layer returns at most `max` rows to one query (the layer's
#' own limit, 1000 or 2000 as a rule): `arcres()` reports that limit and
#' the layer's row count, and a larger layer is read in pages by adding
#' `&resultOffset=` to the address.
#'
#' Errors: a `portal`, `domain`, `layer` or `url` that is not one
#' `http://` or `https://` address (a `domain` may be a bare host name); an
#' `id`, `q` or `key` that is not one non-empty string; a `pid` that is
#' not a table number; a portal that
#' cannot be reached, does not answer in time or answers with an error; an
#' answer that is not what that kind of portal sends; a dataset the portal
#' does not have; a `dest` that is not one file path, or that exists
#' without `over = TRUE`; a `key` `srcls()` does not list.
#'
#' @param portal The address of a CKAN portal, such as
#'   `"https://data.ontario.ca"`.
#' @param id A dataset's identifier on its portal: a CKAN package name or
#'   id, a Socrata four-by-four (`"ijzp-q8t2"`).
#' @param q Search words.
#' @param rows The most datasets `ckanls()` returns, 1 to 1000.
#' @param domain The host of a Socrata portal, such as
#'   `"data.cityofchicago.org"`.
#' @param layer The address of an ArcGIS feature layer, ending in its
#'   number (`.../FeatureServer/0`).
#' @param url The address to download.
#' @param dest The file to write.
#' @param over `TRUE` replaces `dest` when it exists.
#' @param wait The seconds to wait for the file (`getOption("timeout")`
#'   when that is longer).
#' @param key A source's key, as `srcls()` lists it.
#' @return `ckanres()` a data frame with a row per file: `name`, `format`,
#'   `url`, `modified` and `bytes` (`NA` where the portal does not say).
#'   `ckanls()` a data frame with a row per dataset: `id`, `title`, `files`
#'   (their number) and `licence`. `socres()` and `arcres()` a list:
#'   `name`, `url` (the CSV export, or the GeoJSON query), `columns` (their
#'   names), `updated` (a `POSIXct`, or `NA`), and for a layer `rows` and
#'   `max`. `srcget()`, invisibly, a list: `path`, `bytes`, `sha256` and
#'   `url`. `scres()` the address of the table's zipped CSV file.
#'   `srcls()` a data frame: `key`, `publisher`, `what`, `kind` (`"ckan"`,
#'   `"socrata"`, `"arcgis"`, `"statcan"`, `"url"` or `"page"`), `portal`
#'   and `id`. `srcurl()` one address: for a CKAN dataset the first CSV
#'   file's (else the first file's), for a page or a portal its own.
#' @seealso [odtab()] for the tables the project's service mirrors and
#'   checks; [ovrget()] for one oversight report.
#' @examples
#' srcls()[, c("key", "publisher", "kind")]
#' \dontrun{
#' # needs the internet
#' f <- ckanres("https://data.ontario.ca",
#'              "jahn-settlement-data-inmates-ontario")
#' head(f[, c("name", "format", "modified")])
#' s <- socres("data.cityofchicago.org", "dpt3-jri9")
#' s$columns
#' g <- srcget(srcurl("tps-homicides"), tempfile(fileext = ".geojson"))
#' g$sha256
#' scres("35-10-0177-01")
#' }
#' try(srcurl("nowhere"))
#' @export
ckanres <- function(portal, id) {
  fn <- "ckanres"
  portal <- .srweb(portal, "portal", fn)
  id <- .rlstr(id, "id", fn)
  r <- .srckan(
    paste0(portal, "/api/3/action/package_show?id=", .srenc(id)),
    sprintf("the dataset `%s`", id),
    fn
  )
  res <- Filter(is.list, r[["resources"]])
  fld <- function(k) {
    vapply(res, function(x) .srchr(x[[k]]), "", USE.NAMES = FALSE)
  }
  data.frame(
    name = fld("name"),
    format = fld("format"),
    url = fld("url"),
    modified = ifelse(
      is.na(fld("last_modified")),
      fld("created"),
      fld("last_modified")
    ),
    bytes = suppressWarnings(as.numeric(fld("size"))),
    stringsAsFactors = FALSE
  )
}

#' @rdname ckanres
#' @export
ckanls <- function(portal, q, rows = 50L) {
  fn <- "ckanls"
  portal <- .srweb(portal, "portal", fn)
  q <- .rlstr(q, "q", fn)
  if (!.srwhole(rows, 1000)) {
    .ckstop(fn, "`rows` must be a whole number from 1 to 1000")
  }
  r <- .srckan(
    sprintf(
      "%s/api/3/action/package_search?q=%s&rows=%d",
      portal,
      .srenc(q),
      as.integer(rows)
    ),
    "the search",
    fn
  )
  p <- Filter(is.list, r[["results"]])
  data.frame(
    id = vapply(p, function(x) .srchr(x[["name"]]), ""),
    title = vapply(p, function(x) .srchr(x[["title"]]), ""),
    files = vapply(p, function(x) length(x[["resources"]]), 0L),
    licence = vapply(p, function(x) .srchr(x[["license_title"]]), ""),
    stringsAsFactors = FALSE
  )
}

#' @rdname ckanres
#' @export
socres <- function(domain, id) {
  fn <- "socres"
  base <- .srweb(.srhost(domain), "domain", fn)
  id <- .rlstr(id, "id", fn)
  if (!grepl("^[a-z0-9]{4}-[a-z0-9]{4}$", id)) {
    .ckstop(fn, "`id` must be a Socrata identifier, such as \"ijzp-q8t2\"")
  }
  m <- .srjson(
    paste0(base, "/api/views/", id, ".json"),
    sprintf("the dataset `%s`", id),
    fn
  )
  if (!is.list(m) || !identical(.srchr(m[["id"]]), id)) {
    .ckstop(fn, "%s did not describe the dataset `%s`", base, id)
  }
  cols <- Filter(is.list, m[["columns"]])
  up <- suppressWarnings(as.numeric(m[["rowsUpdatedAt"]]))
  list(
    name = .srchr(m[["name"]]),
    url = sprintf("%s/api/views/%s/rows.csv?accessType=DOWNLOAD", base, id),
    columns = vapply(cols, function(x) .srchr(x[["name"]]), ""),
    updated = .srtime(if (length(up) == 1L) up else NA_real_)
  )
}

#' @rdname ckanres
#' @export
arcres <- function(layer) {
  fn <- "arcres"
  layer <- .srweb(layer, "layer", fn)
  if (!grepl("/(FeatureServer|MapServer)/[0-9]+$", layer)) {
    .ckstop(
      fn,
      "`layer` must end in a layer number, as .../FeatureServer/0"
    )
  }
  m <- .srjson(paste0(layer, "?f=json"), "the layer", fn)
  if (!is.list(m) || !is.null(m[["error"]]) || is.null(m[["name"]])) {
    .ckstop(fn, "%s did not describe a layer", layer)
  }
  n <- .srjson(
    sprintf("%s/query?where=1%%3D1&returnCountOnly=true&f=json", layer),
    "the layer's row count",
    fn
  )
  num <- function(v) {
    v <- suppressWarnings(as.numeric(v))
    if (length(v) == 1L) v else NA_real_
  }
  ed <- m[["editingInfo"]]
  list(
    name = .srchr(m[["name"]]),
    url = sprintf("%s/query?where=1%%3D1&outFields=*&f=geojson", layer),
    columns = vapply(
      Filter(is.list, m[["fields"]]),
      function(x) .srchr(x[["name"]]),
      ""
    ),
    updated = .srtime(
      num(if (is.list(ed)) ed[["lastEditDate"]] else NULL) / 1000
    ),
    rows = num(if (is.list(n)) n[["count"]] else NULL),
    max = num(m[["maxRecordCount"]])
  )
}

#' @rdname ckanres
#' @export
srcget <- function(url, dest, over = FALSE, wait = 600) {
  fn <- "srcget"
  url <- .srweb(url, "url", fn, file = TRUE)
  over <- .ckflag(over, "over", fn)
  .srdest(dest, wait, over, fn)
  tmp <- paste0(dest, ".part")
  .odget(url, tmp, wait, "the file", fn)
  if (!file.rename(tmp, dest)) {
    unlink(tmp)
    .ckstop(fn, "%s could not be written", dest)
  }
  invisible(list(
    path = dest,
    bytes = file.size(dest),
    sha256 = .odhash(dest, fn),
    url = url
  ))
}

#' @rdname ckanres
#' @export
srcls <- function() {
  r <- .srreg()
  as.data.frame(
    stats::setNames(
      lapply(seq_len(ncol(r)), function(j) unname(r[, j])),
      colnames(r)
    ),
    stringsAsFactors = FALSE
  )
}

#' @rdname ckanres
#' @export
srcurl <- function(key) {
  fn <- "srcurl"
  key <- .rlstr(key, "key", fn)
  r <- .srreg()
  i <- match(key, r[, "key"])
  if (is.na(i)) {
    .ckstop(fn, "`%s` is not a key of srcls()", key)
  }
  s <- r[i, ]
  if (!nzchar(s[["id"]])) {
    return(s[["portal"]])
  }
  switch(
    s[["kind"]],
    statcan = .clas(fn, scres(s[["id"]])),
    ckan = {
      f <- .clas(fn, ckanres(s[["portal"]], s[["id"]]))
      if (!nrow(f)) {
        .ckstop(fn, "%s lists no file for `%s`", s[["portal"]], s[["id"]])
      }
      csv <- which(toupper(f[["format"]]) == "CSV")
      f[["url"]][if (length(csv)) csv[1L] else 1L]
    },
    socrata = .clas(fn, socres(s[["portal"]], s[["id"]]))[["url"]],
    arcgis = .clas(
      fn,
      arcres(sprintf("%s/%s/FeatureServer/0", s[["portal"]], s[["id"]]))
    )[["url"]]
  )
}

#' @rdname ckanres
#' @param pid The number of a Statistics Canada table: eight digits, with
#'   or without the hyphens of `"35-10-0177-01"` (the last two digits, the
#'   view, are not part of it).
#' @param lang `"en"` (the default) or `"fr"`: the language of the file.
#' @export
scres <- function(pid, lang = c("en", "fr")) {
  fn <- "scres"
  lang <- .ckmarg("lang", fn)
  if (is.numeric(pid) && length(pid) == 1L && is.finite(pid)) {
    pid <- format(pid, scientific = FALSE)
  }
  pid <- gsub("-", "", .rlstr(pid, "pid", fn), fixed = TRUE)
  if (!grepl("^[0-9]{8}([0-9]{2})?$", pid)) {
    .ckstop(
      fn,
      "`pid` must be a table number of eight digits, such as \"35100177\""
    )
  }
  pid <- substr(pid, 1L, 8L)
  r <- .srjson(
    sprintf(
      "https://www150.statcan.gc.ca/t1/wds/rest/getFullTableDownloadCSV/%s/%s",
      pid,
      lang
    ),
    sprintf("the table `%s`", pid),
    fn
  )
  u <- .scurl(r)
  if (is.na(u)) {
    .ckstop(fn, "Statistics Canada did not find the table `%s`", pid)
  }
  u
}

#' The tables of the data service about Indigenous peoples
#'
#' Lists the tables of the project's open data service drawn from
#' Statistics Canada's tables on Indigenous peoples and from the 2021
#' Census Indigenous Population Profile, with their descriptions, sizes
#' and licences, optionally those whose description matches a pattern.
#' [odtab()] loads one by its key.
#'
#' Only openly licensed tables of aggregates are mirrored. The data of
#' the First Nations Information Governance Centre (the First Nations
#' Regional Health Survey and its companions) are governed by the OCAP
#' principles of ownership, control, access and possession, and are not;
#' they are asked for from the Centre.
#'
#' Errors: `pat` not one non-empty string or not a valid regular
#' expression; those of [odls()].
#'
#' @param pat A regular expression matched against the descriptions
#'   without regard to case, or `NULL` (all tables).
#' @inheritParams odtab
#' @return A data frame with a row per table: `key`, `rows`,
#'   `description`, `licence` and `source`.
#' @seealso [odtab()] and [odls()]; [indiex()] for the tables by subject.
#' @examples
#' \dontrun{
#' # needs the internet
#' t <- indtabs("mental health")
#' t[, c("key", "rows")]
#' }
#' try(indtabs(pat = ""))
#' @export
indtabs <- function(
  pat = NULL,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
) {
  fn <- "indtabs"
  if (!is.null(pat)) {
    pat <- .rlstr(pat, "pat", fn)
    if (inherits(try(grepl(pat, ""), silent = TRUE), "try-error")) {
      .ckstop(fn, "`pat` is not a valid regular expression")
    }
  }
  base <- .odbase(base, fn)
  ds <- Filter(
    function(d) {
      is.list(d) &&
        is.character(d[["key"]]) &&
        grepl("^(statcan_indigenous|census2021_ipp)/", d[["key"]])
    },
    .odman(base, fn)[["datasets"]]
  )
  meta <- function(k) {
    vapply(
      ds,
      function(d) {
        m <- d[["meta"]]
        v <- if (is.list(m)) .srchr(m[[k]]) else NA_character_
        if (is.na(v)) "" else v
      },
      "",
      USE.NAMES = FALSE
    )
  }
  out <- data.frame(
    key = vapply(ds, function(d) d[["key"]], ""),
    rows = .odfld(ds, "rows", NA_real_),
    description = meta("description"),
    licence = meta("license"),
    source = meta("source"),
    stringsAsFactors = FALSE
  )
  if (!is.null(pat)) {
    out <- out[
      grepl(pat, out[["description"]], ignore.case = TRUE),
      ,
      drop = FALSE
    ]
  }
  rownames(out) <- NULL
  out
}

#' Statistics Canada tables on Indigenous peoples, by subject
#'
#' The tables of the project's open data service that hold what is most
#' often asked for about First Nations people, Metis and Inuit: health
#' indicators, perceived mental health, wellbeing, suicidal thoughts, the
#' census counts of population and of households, life expectancy and
#' infant mortality. `indiex()` finds the tables of one subject among
#' those [indtabs()] lists, and with `load = TRUE` loads the first.
#'
#' The package holds no data: the tables are Statistics Canada's, mirrored
#' by the service under the Statistics Canada Open Licence and loaded and
#' checked by [odtab()]. A subject is matched by words in a table's
#' description, so a table the service adds later is found without a new
#' release; a subject with no table there gives a data frame of no rows.
#'
#' Errors: `name` not one of the subjects; `load` not `TRUE` or `FALSE`;
#' with `load`, a subject with no table; those of [indtabs()] and
#' [odtab()].
#'
#' @param name The subject: `"health"`, `"mentalhealth"`, `"wellbeing"`,
#'   `"suicide"`, `"censuspop"`, `"censushh"`, `"lifeexp"` or
#'   `"infantmort"`.
#' @param load `TRUE` loads the first table found, with [odtab()].
#' @inheritParams odtab
#' @return The rows of [indtabs()] for the subject; with `load = TRUE`,
#'   the first of those tables as a data frame.
#' @seealso [indtabs()], [odtab()]; [sctidy()] and [scests()] to read a
#'   Statistics Canada table once loaded.
#' @examples
#' \dontrun{
#' # needs the internet
#' indiex("mentalhealth")[, c("key", "rows")]
#' x <- indiex("lifeexp", load = TRUE)
#' }
#' try(indiex("weather"))
#' @export
indiex <- function(
  name = c(
    "health",
    "mentalhealth",
    "wellbeing",
    "suicide",
    "censuspop",
    "censushh",
    "lifeexp",
    "infantmort"
  ),
  load = FALSE,
  base = getOption("crimlab.data", "https://data.crimlab.dev")
) {
  fn <- "indiex"
  name <- .ckmarg("name", fn)
  load <- .ckflag(load, "load", fn)
  pat <- c(
    health = "health indicator|chronic condition|health characteristic",
    mentalhealth = "mental health",
    wellbeing = "well-?being|life satisfaction",
    suicide = "suicid",
    censuspop = "population|identity",
    censushh = "household|dwelling",
    lifeexp = "life expectancy",
    infantmort = "infant mortal"
  )[[name]]
  t <- .clas(fn, indtabs(pat, base = base))
  if (name %in% c("censuspop", "censushh")) {
    t <- t[startsWith(t[["key"]], "census2021_ipp/"), , drop = FALSE]
    rownames(t) <- NULL
  }
  if (!load) {
    return(t)
  }
  if (!nrow(t)) {
    .ckstop(fn, "the service holds no table on `%s`", name)
  }
  .clas(fn, odtab(t[["key"]][1L], base = base))
}

#' Fetch one police oversight report
#'
#' Downloads a director's report of Ontario's Special Investigations
#' Unit by its report number, as the page the Unit serves, to a file.
#' [ovrtxt()] then gives its text and [ovrpar()] codes it into fields.
#'
#' The report number is the `drid` in the address of a report's page. The
#' Unit publishes each report in English and in French.
#'
#' Errors: `drid` not one positive whole number; `lang` not `"en"` or
#' `"fr"`; those of [srcget()].
#'
#' @param drid The report number: one positive whole number.
#' @param lang `"en"` (the default) or `"fr"`.
#' @inheritParams ckanres
#' @return As [srcget()]: invisibly, a list of `path`, `bytes`, `sha256`
#'   and `url`.
#' @seealso [ovrtxt()], [ovrpar()]; [srcls()] for other sources.
#' @examples
#' \dontrun{
#' # needs the internet
#' f <- ovrget(1234, tempfile(fileext = ".html"))
#' substr(ovrtxt(readLines(f$path, warn = FALSE)), 1, 200)
#' }
#' try(ovrget(0, tempfile()))
#' @export
ovrget <- function(drid, dest, lang = c("en", "fr"), over = FALSE, wait = 120) {
  fn <- "ovrget"
  lang <- .ckmarg("lang", fn)
  if (!.srwhole(drid, .Machine$integer.max)) {
    .ckstop(
      fn,
      paste(
        "`drid` must be one positive report number",
        "(the drid of the report's address)"
      )
    )
  }
  .clas(
    fn,
    srcget(
      sprintf(
        "https://www.siu.on.ca/%s/directors_report_details.php?drid=%d",
        lang,
        as.integer(drid)
      ),
      dest,
      over = over,
      wait = wait
    )
  )
}

#' The sources srcls() lists: key, publisher, what, kind, portal, id. Each
#' was resolved against its portal when it was added; the test of the
#' registry does so again when run with the network. A source with no id
#' is a portal to search (ckanls()) or a page to read.
#' @noRd
.srnm <- list(
  tps = "https://services.arcgis.com/S9th0jAJ7bqgIRjw/arcgis/rest/services",
  tor = "https://ckan0.cf.opendata.inter.prod-toronto.ca",
  ont = "https://data.ontario.ca",
  sc = "https://www150.statcan.gc.ca",
  es = "https://ec.europa.eu/eurostat/api/dissemination/sdmx/2.1/data/",
  tp = "Toronto Police Service",
  sg = "Ontario Ministry of the Solicitor General",
  st = "Statistics Canada",
  ny = "New York City Police Department",
  ch = "Chicago Police Department"
)

#' The sources of srcls(): six fields a source (.srnm holds the names
#' that repeat).
#' @noRd
.srtab <- with(.srnm, {
  x <- c(
    # Canada: police
    "tps-mci",
    tp,
    "Major crime indicators",
    "arcgis",
    tps,
    "Major_Crime_Indicators_Open_Data",
    "tps-homicides",
    tp,
    "Homicides",
    "arcgis",
    tps,
    "Homicides_Open_Data_ASR_RC_TBL_002",
    "tps-shootings",
    tp,
    "Shootings and firearm discharges",
    "arcgis",
    tps,
    "Shooting_and_Firearm_Discharges_Open_Data",
    "tps-rates",
    tp,
    "Neighbourhood crime rates",
    "arcgis",
    tps,
    "Neighbourhood_Crime_Rates_Open_Data",
    "tps-uof",
    tp,
    "Race-based data: use of force",
    "arcgis",
    tps,
    "RBDC_UOF_TBL_001",
    "tps-arrests",
    tp,
    "Race-based data: arrests and strip searches",
    "arcgis",
    tps,
    "RBDC_ARR_TBL_001",
    "tps-mha",
    tp,
    "Mental Health Act apprehensions",
    "arcgis",
    tps,
    "Mental_Health_Act_Apprehensions_Open_Data",
    "toronto-asr-crimes",
    "City of Toronto",
    "Police annual statistical report: reported crimes",
    "ckan",
    tor,
    "police-annual-statistical-report-reported-crimes",
    "toronto-asr-arrests",
    "City of Toronto",
    "Police annual statistical report: arrested and charged persons",
    "ckan",
    tor,
    "police-annual-statistical-report-arrested-and-charged-persons",
    "montreal-crimes",
    "Ville de Montreal",
    "Criminal acts (actes criminels)",
    "ckan",
    "https://donnees.montreal.ca",
    "actes-criminels",
    "siu-reports",
    "Special Investigations Unit of Ontario",
    "Director's reports, one page each (ovrget())",
    "page",
    "https://www.siu.on.ca/en/directors_reports.php",
    "",
    "vpd-crime",
    "Vancouver Police Department",
    "Crime data, behind a form that accepts its terms",
    "page",
    "https://geodash.vpd.ca/opendata/",
    "",
    # Canada: corrections and courts
    "ontario-jahn",
    sg,
    "Jahn settlement: data on people in custody",
    "ckan",
    ont,
    "jahn-settlement-data-inmates-ontario",
    "ontario-uof",
    sg,
    "Use of force in correctional institutions",
    "ckan",
    ont,
    "use-of-force-in-correctional-institutions",
    "ontario-supervision",
    sg,
    "Profile of people under correctional supervision",
    "ckan",
    ont,
    "people-under-correctional-supervision",
    # Canada: Statistics Canada tables
    "statcan-crime",
    st,
    "Incident-based crime statistics, by detailed violations",
    "statcan",
    sc,
    "35100177",
    "statcan-csi",
    st,
    "Crime severity index and weighted clearance rates",
    "statcan",
    sc,
    "35100026",
    "statcan-homicide",
    st,
    "Homicide victims: number and rate",
    "statcan",
    sc,
    "35100068",
    "statcan-hate",
    st,
    "Police-reported hate crime",
    "statcan",
    sc,
    "35100191",
    "statcan-courts",
    st,
    "Adult criminal courts: cases and charges by type of decision",
    "statcan",
    sc,
    "35100027",
    "statcan-admissions",
    st,
    "Adult admissions to correctional services",
    "statcan",
    sc,
    "35100014",
    "statcan-custody",
    st,
    "Average counts of adults in correctional programs",
    "statcan",
    sc,
    "35100154",
    # United States
    "nypd-complaints",
    ny,
    "Complaints, historic",
    "socrata",
    "data.cityofnewyork.us",
    "qgea-i56i",
    "nypd-complaints-ytd",
    ny,
    "Complaints, year to date",
    "socrata",
    "data.cityofnewyork.us",
    "5uac-w243",
    "nypd-arrests",
    ny,
    "Arrests, historic",
    "socrata",
    "data.cityofnewyork.us",
    "8h9b-rp9u",
    "nypd-arrests-ytd",
    ny,
    "Arrests, year to date",
    "socrata",
    "data.cityofnewyork.us",
    "uip8-fykc",
    "nypd-hate",
    ny,
    "Hate crimes",
    "socrata",
    "data.cityofnewyork.us",
    "bqiq-cu78",
    "cpd-crimes",
    ch,
    "Crimes, 2001 to present",
    "socrata",
    "data.cityofchicago.org",
    "ijzp-q8t2",
    "cpd-arrests",
    ch,
    "Arrests",
    "socrata",
    "data.cityofchicago.org",
    "dpt3-jri9",
    "lapd-crimes",
    "Los Angeles Police Department",
    "Crime data, 2020 to 2024",
    "socrata",
    "data.lacity.org",
    "2nrs-mtv8",
    "fbi-cde",
    "Federal Bureau of Investigation",
    "Crime Data Explorer (downloads and an interface that needs a key)",
    "page",
    "https://cde.ucr.cjis.gov/LATEST/webapp/",
    "",
    "bjs-ncvs",
    "Bureau of Justice Statistics",
    "National Crime Victimization Survey",
    "page",
    "https://bjs.ojp.gov/data-collection/ncvs",
    "",
    # United Kingdom, Europe, Australia
    "uk-police",
    "Home Office (data.police.uk)",
    "Street-level crime, outcomes and stop and search: the latest archive",
    "url",
    "https://data.police.uk/data/archive/latest.zip",
    "",
    "eurostat-crime",
    "Eurostat",
    "Police-recorded offences by offence category (crim_off_cat)",
    "url",
    paste0(es, "crim_off_cat?format=TSV&compressed=false"),
    "",
    "nsw-bocsar",
    "NSW Bureau of Crime Statistics and Research",
    "Open datasets",
    "page",
    "https://bocsar.nsw.gov.au/statistics-dashboards/open-datasets.html",
    "",
    # portals to search with ckanls()
    "portal-canada",
    "Government of Canada",
    "Open Government portal",
    "ckan",
    "https://open.canada.ca/data",
    "",
    "portal-ontario",
    "Government of Ontario",
    "Ontario data catalogue",
    "ckan",
    ont,
    "",
    "portal-toronto",
    "City of Toronto",
    "Open data portal",
    "ckan",
    tor,
    "",
    "portal-montreal",
    "Ville de Montreal",
    "Donnees ouvertes",
    "ckan",
    "https://donnees.montreal.ca",
    "",
    "portal-uk",
    "Government of the United Kingdom",
    "data.gov.uk",
    "ckan",
    "https://data.gov.uk",
    "",
    "portal-australia",
    "Australian Government",
    "data.gov.au",
    "ckan",
    "https://data.gov.au/data",
    ""
  )
  matrix(
    x,
    ncol = 6L,
    byrow = TRUE,
    dimnames = list(
      NULL,
      c("key", "publisher", "what", "kind", "portal", "id")
    )
  )
})

#' The table of sources.
#' @noRd
.srreg <- function() .srtab

#' One http(s):// address (or, with file, a file:// one) without trailing
#' slashes, or the refusal of fn.
#' @noRd
.srweb <- function(x, what, fn, file = FALSE) {
  re <- if (file) {
    "^(https?|file)://[^[:space:]]+$"
  } else {
    "^https?://[^[:space:]]+$"
  }
  if (
    !is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !validEnc(x) ||
      !grepl(re, x)
  ) {
    .ckstop(
      fn,
      "`%s` must be one address starting with http:// or https://",
      what
    )
  }
  sub("/+$", "", x)
}

#' Text as part of a query: everything but letters, digits, - . _ ~ as
#' percent escapes.
#' @noRd
.srenc <- function(x) utils::URLencode(x, reserved = TRUE)

#' One string from a parsed JSON value, NA when it is not one.
#' @noRd
.srchr <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    x
  } else if (is.numeric(x) && length(x) == 1L) {
    format(x, scientific = FALSE)
  } else {
    NA_character_
  }
}

#' Seconds since 1970 as a time in UTC, NA kept.
#' @noRd
.srtime <- function(s) {
  as.POSIXct(s, origin = "1970-01-01", tz = "UTC")
}

#' The JSON a portal answers at url, as lists; a failure refuses in fn's
#' name, saying what was asked for.
#' @noRd
.srjson <- function(url, what, fn) {
  tmp <- tempfile("clsrc")
  on.exit(unlink(tmp), add = TRUE)
  .odget(url, tmp, 60, what, fn)
  txt <- readBin(tmp, "raw", file.size(tmp))
  tryCatch(
    .odjson(txt),
    error = function(e) {
      .ckstop(
        fn,
        "%s answered %s with something that is not JSON",
        sub("^(https?://[^/]*).*$", "\\1", url),
        what
      )
    }
  )
}

#' The `result` of a CKAN action at url; a portal that answers without
#' success refuses in fn's name.
#' @noRd
.srckan <- function(url, what, fn) {
  r <- .srjson(url, what, fn)
  if (!is.list(r) || !isTRUE(r[["success"]]) || !is.list(r[["result"]])) {
    .ckstop(
      fn,
      "%s did not find %s",
      sub("^(https?://[^/]*).*$", "\\1", url),
      what
    )
  }
  r[["result"]]
}

#' dest and wait of srcget() checked; a file in the way is refused.
#' @noRd
.srdest <- function(dest, wait, over, fn) {
  if (
    !is.character(dest) || length(dest) != 1L || is.na(dest) || !nzchar(dest)
  ) {
    .ckstop(fn, "`dest` must be one file path")
  }
  if (
    !is.numeric(wait) || length(wait) != 1L || !is.finite(wait) || wait <= 0
  ) {
    .ckstop(fn, "`wait` must be one positive number of seconds")
  }
  if (file.exists(dest) && !over) {
    .ckstop(fn, "%s exists; `over = TRUE` replaces it", dest)
  }
  invisible(NULL)
}


#' One whole number from 1 to `hi`?
#' @noRd
.srwhole <- function(v, hi) {
  is.numeric(v) &&
    length(v) == 1L &&
    is.finite(v) &&
    v >= 1 &&
    v <= hi &&
    v == floor(v)
}

#' A host name given without its scheme takes https.
#' @noRd
.srhost <- function(domain) {
  bare <- is.character(domain) &&
    length(domain) == 1L &&
    !is.na(domain) &&
    !grepl("^https?://", domain)
  if (bare) paste0("https://", domain) else domain
}

#' The address in Statistics Canada's answer `r`, or NA when it found no
#' table.
#' @noRd
.scurl <- function(r) {
  if (!is.list(r) || !identical(.srchr(r[["status"]]), "SUCCESS")) {
    return(NA_character_)
  }
  u <- .srchr(r[["object"]])
  if (is.na(u) || !grepl("^https://", u)) NA_character_ else u
}
