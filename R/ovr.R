#' Text of a police oversight report
#'
#' Turns the page of a report published by a police oversight body (a
#' civilian investigation unit, an ombudsman, an inspector) into the plain
#' text its other tools read: markup removed, each paragraph, heading,
#' list item and table row on a line of its own, so that section headings
#' stand alone, and character entities decoded. Use it to read or search
#' saved reports without a browser, or before the rules that code them.
#'
#' Scripts, styles and comments are dropped. The start and end of a block
#' element (a paragraph, a heading, a division, a list item, a table row,
#' a line break and the like) end a line; every other tag is a space, and
#' the page's own line breaks are spaces, as a browser shows them. Named
#' entities common on report pages and every numeric entity are decoded
#' (curly quotes and dashes to their plain forms). Spaces are collapsed,
#' lines trimmed, and at most one empty line is kept between paragraphs.
#' Plain text (with `html = FALSE`, or with no markup when `html = NA`) is
#' only tidied that way, its line breaks kept. Nothing is fetched and
#' nothing is written; the pass is linear in the length of the page.
#'
#' @param x The pages or texts: a character vector. `NA` stays `NA`.
#' @param html `NA` (the default) to read each element as HTML when it
#'   holds markup (a `<` followed by a letter, `/` or `!`), `TRUE` or
#'   `FALSE` to say.
#' @return A character vector the length of `x`, in UTF-8.
#' @examples
#' ovrtxt("<h2>The Team</h2><p>Number of investigators assigned: 3</p>")
#' ovrtxt("<p>WO\n#1 Interviewed</p><p>WO&nbsp;#2 Interviewed</p>")
#' ovrtxt("a  b\r\nc", html = FALSE)
#' @export
ovrtxt <- function(x, html = NA) {
  f <- "ovrtxt"
  x <- .ovrin(x, "x", f)
  .Call(C_clab_r_ovrtxt, x, .ovrhtml(html, f))
}

#' Dates written out in reports, as ISO dates
#'
#' Finds the first date written in each text, in English or French, and
#' gives it as `YYYY-MM-DD`, the form that sorts and that [as.Date()]
#' reads. `"January 5, 2023"`, `"Jan. 5th 2023"`, `"5 January 2023"`,
#' `"5th of January, 2023"`, `"5 janvier 2023"` and `"2023-01-05"` all
#' give `"2023-01-05"`.
#'
#' Month names are English or French, whole or abbreviated (with or
#' without a full stop), in any case, French accents optional; the day
#' may carry an ordinal (`"1st"`, `"1er"`). A day that does not exist
#' (`"February 30, 2023"`) is passed over, and the search goes on to the
#' next date in the text. The reading never depends on the session's
#' locale.
#'
#' @param x The texts: a character vector.
#' @return A character vector the length of `x`: the first date of each
#'   as `YYYY-MM-DD`, `NA` where there is none (or `x` is `NA`).
#' @examples
#' ovrdte(c("January 5, 2023", "le 1er février 2021", "no date"))
#' ovrdte("On February 30, 2023 or rather March 2, 2023, the SIU ...")
#' as.Date(ovrdte("Date: April 28, 2023"))
#' @export
ovrdte <- function(x) {
  .Call(C_clab_r_ovrdte, .ovrin(x, "x", "ovrdte"))
}

#' Field schemas for coding police oversight reports
#'
#' The fields a preset reads from each report, and the rule, cues and
#' section of each: a data frame that [ovrpar()] takes as `schema`. Copy
#' it and change rows, or write one of the same columns, to code the
#' reports of any oversight body; nothing in the rules is tied to one
#' agency.
#'
#' A field is read by one rule. `cue` holds the phrases the rule looks for,
#' alternatives separated by `"|"`, matched as written (case counts) and
#' as whole words at their ends. With `sect` (headings separated by
#' `"|"`), the rule reads only the section under the last line that is
#' one of those headings (a page that opens with a table of contents
#' repeats each heading), up to the next line that looks like a heading
#' (a short line starting with a capital, with no digit, not ending in a
#' stop or a comma); a report without the heading leaves the field
#' missing. The rules:
#'
#' * `after`: the text after a cue (past spaces and a colon), up to a
#'   comma, a semicolon, a full stop or the end of the line.
#' * `num`: the first number (digits, or a number word in English or
#'   French up to twenty) within 40 characters after a cue on its line.
#' * `date`: the first date ([ovrdte()]) in the sentence of a cue, or
#'   with no cue the first date in the section.
#' * `tag`: the highest number of a numbered tag (`"WO #3"`, `"WO 3"`,
#'   `"AT no 3"`): the count of the people the tags number.
#' * `sent`: the sentence that holds the first cue.
#' * `prev`: the line before the last line that is a cue (a signature:
#'   the name above a title).
#' * `lines`: every line that starts with a cue, distinct, joined by
#'   `"; "`.
#' * `name`: a name of up to four capitalised words that ends in a cue
#'   (`"Barrie Police Service"` for the cue `"Police Service"`), without
#'   a leading "The".
#' * `flag`: `FALSE` when a phrase of `alt` occurs, else `TRUE` when a
#'   cue does.
#' * `before`: the number just before a cue (`"34"` in `"34-year-old"`).
#' * `word`: the word just after a cue (`"man"` in `"34-year-old man"`).
#'
#' The preset `"siu"` reads the director's reports of Ontario's Special
#' Investigations Unit as published in English since 2020: the case
#' number, the police service, the dates of the incident, of the
#' notification and of the decision, the investigators assigned, the
#' witness officials, civilian witnesses and subject officials, the age
#' and sex of the affected person, whether the director found grounds
#' for charges, the director's name, the municipality, the injuries and
#' the provisions considered.
#'
#' @param preset The name of a preset: `"siu"`.
#' @return A data frame with one row per field: `name`, `rule`, `cue`,
#'   `alt`, `sect` and `desc` (what the field records).
#' @examples
#' s <- ovrsch()
#' s[, c("name", "rule")]
#' # a field of one's own: the date a complaint came in
#' s2 <- rbind(s, data.frame(name = "date_received", rule = "date",
#'                           cue = "received the complaint", alt = "",
#'                           sect = "",
#'                           desc = "the date the complaint came in"))
#' tail(s2[, c("name", "rule", "cue")], 2)
#' @export
ovrsch <- function(preset = "siu") {
  f <- "ovrsch"
  preset <- .ckarg(preset, "siu", "preset", f)
  as.data.frame(.Call(C_clab_r_ovrsch, preset), stringsAsFactors = FALSE)
}

#' Code police oversight reports into fields
#'
#' Reads the fields of a schema from each report: the dates, the counts of
#' officials and witnesses, the service, the outcome and whatever else
#' the schema names, each by its rule (see [ovrsch()]). The rules are
#' deterministic and offline, so the same report always gives the same
#' row, and a folder of saved reports becomes a data set whose every
#' value can be traced to the words that gave it.
#'
#' Each report is first made text as [ovrtxt()] makes it, then each field
#' is read by its rule; a field the report does not state is `NA`. A
#' field of the rules `num`, `tag` and `before` is a whole number, one of
#' the rule `flag` `TRUE` or `FALSE`, one of the rule `date` a
#' `YYYY-MM-DD` date (as text), and the others text. Nothing is fetched.
#'
#' @inheritParams ovrtxt
#' @param x The reports: their pages (HTML) or texts, a character vector.
#' @param schema The fields: a data frame with the character columns
#'   `name`, `rule`, `cue`, `alt` and `sect`, as [ovrsch()] gives them.
#' @return A data frame with one row per report and one column per field,
#'   in the order of `schema`; a row of `NA` for a report that is `NA`.
#' @examples
#' r <- paste(
#'   "<h1>Director's Report - Case # 23-OFD-001</h1>",
#'   "<p>On January 6, 2023, the Barrie Police Service contacted the SIU.",
#'   "In the night, a 34-year-old man was injured.</p>",
#'   "<p>Number of SIU Investigators assigned: 3</p>",
#'   "<p>WO #1 Interviewed</p><p>WO #2 Interviewed</p>",
#'   "<p>Joseph Martino</p><p>Director</p>"
#' )
#' p <- ovrpar(r)
#' p[, c("case_number", "police_service", "date_siu_notified_iso",
#'       "siu_investigators", "number_of_witness_officials",
#'       "age_affected", "directors_name")]
#' @export
ovrpar <- function(x, schema = ovrsch(), html = NA) {
  f <- "ovrpar"
  x <- .ovrin(x, "x", f)
  html <- .ovrhtml(html, f)
  s <- .ovrschk(schema, f)
  r <- .Call(C_clab_r_ovrpar, x, html, match(s$rule, .ovrrules) - 1L,
             s$cue, s$alt, s$sect)
  r <- lapply(seq_along(r), function(j) {
    switch(s$rule[j],
           num = , tag = , before = as.integer(r[[j]]),
           flag = as.logical(r[[j]]),
           r[[j]])
  })
  as.data.frame(stats::setNames(r, s$name), stringsAsFactors = FALSE,
                check.names = FALSE)
}

#' The number of subject officials in an oversight report
#'
#' Counts the subject officials of each report, the officers whose
#' conduct the oversight body investigated, and says which rule gave the
#' count and on what words, so that the number can be checked against the
#' report. Use it where a report has no roster [ovrpar()] can read, or to
#' check a count coded by hand.
#'
#' The rules run in order and the first that applies gives the count:
#' 1. the roster under a heading of `head`: the highest numbered tag
#'    there (`"SO #2"`);
#' 2. the highest numbered tag in the whole report, trusted only when the
#'    tag numbered 1 is there too;
#' 3. a number before a plural of `term` (`"two subject officials"`,
#'    `"3 subject officers"`);
#' 4. a mention of `term` in the singular, not after "no" or "not a": 1;
#' 5. a statement that there was none (`"no subject official"`), or
#'    numbered witness tags (`wtag`) with no subject tag at all: 0, a
#'    real answer;
#' 6. otherwise unresolved: `NA`.
#'
#' Tags and headings are matched as written; `term` in any case. The
#' defaults are the words of Canadian oversight bodies; give a body's own
#' (`tag = "AI"`, `term = "agent impliqué"` for a French report).
#'
#' @inheritParams ovrtxt
#' @param x The reports: pages (HTML) or texts, a character vector.
#' @param tag The tag of a subject official in the report: one text.
#' @param term The words for a subject official, singular: text.
#' @param wtag The tag of a witness official: one text.
#' @param head The headings of the roster of subject officials: text.
#' @return A data frame with one row per report: `count` (whole number,
#'   `NA` when unresolved), `rule` (1 to 6, as above) and `reason` (the
#'   rule and the words it read).
#' @examples
#' ovrsoc(c("Subject Officials\nSO #1 Interviewed\nSO #2 Declined",
#'          "Two subject officials were designated.",
#'          "Witness Officials\nWO #1\nWO #2",
#'          "No officers are named here."))
#' @export
ovrsoc <- function(x, tag = "SO", term = c("subject official",
                                           "subject officer"),
                   wtag = "WO", head = c("Subject Officials",
                                         "Subject Officers"),
                   html = NA) {
  f <- "ovrsoc"
  x <- .ovrin(x, "x", f)
  a <- list(tag = tag, term = term, wtag = wtag, head = head)
  for (nm in names(a)) {
    a[[nm]] <- .ovrwords(a[[nm]], nm, f, one = nm %in% c("tag", "wtag"))
  }
  html <- .ovrhtml(html, f)
  r <- .Call(C_clab_r_ovrsoc, x, html, a$tag, a$term, a$wtag, a$head)
  as.data.frame(r, stringsAsFactors = FALSE)
}

#' Keyword indicators in oversight reports
#'
#' Marks each report that mentions any word of each keyword group: the
#' signals of mental health and of race that analyses of oversight
#' reports tally, or any groups of one's own. A mention is a signal to
#' read the report, not a finding: a report may name a condition only to
#' rule it out.
#'
#' Words are matched whole and in any case (ASCII letters; others as
#' written); a word ending in `*` matches any word it starts
#' (`"schizophren*"` finds "schizophrenia"). The preset `"mhr"` has the
#' groups `mental_health` (mental health or illness, psychiatric,
#' schizophrenia, bipolar, depression, psychosis, suicidal, self-harm,
#' delusions, crisis intervention, "EDP", "emotionally disturbed",
#' "MCIT") and `race` (Black, African, Indigenous, First Nations, Metis,
#' Inuit, racialized, racial, racism, anti-Black, ethnic, South Asian).
#'
#' @inheritParams ovrtxt
#' @param x The reports: pages (HTML) or texts, a character vector.
#' @param vocab The keyword groups: the name of a preset (`"mhr"`), a
#'   named list of character vectors (one group each), or a character
#'   vector (each word a group of its own).
#' @return A data frame with one row per report: a logical column per
#'   group (`TRUE` when the report holds a word of it; `NA` for a report
#'   that is `NA`) and `terms`, the words found, in the order of `vocab`,
#'   separated by `"; "`.
#' @examples
#' x <- c("The man had a history of schizophrenia; a crisis team came.",
#'        "An Indigenous woman was arrested.",
#'        "Nothing of either kind.")
#' k <- ovrkwd(x)
#' k
#' colSums(k[c("mental_health", "race")])
#' ovrkwd(x, list(force = c("taser", "conducted energy weapon", "baton")))
#' @export
ovrkwd <- function(x, vocab = "mhr", html = NA) {
  f <- "ovrkwd"
  x <- .ovrin(x, "x", f)
  v <- .ovrvoc(vocab, f)
  html <- .ovrhtml(html, f)
  m <- .Call(C_clab_r_ovrkwd, x, html, v$term)
  hit <- m > 0L
  gs <- unique(v$group)
  out <- lapply(gs, function(g) {
    h <- hit[, v$group == g, drop = FALSE]
    rowSums(h) > 0L
  })
  out <- as.data.frame(stats::setNames(out, gs), check.names = FALSE)
  out$terms <- vapply(seq_along(x), function(i) {
    if (is.na(x[i])) {
      return(NA_character_)
    }
    paste(v$term[hit[i, ]], collapse = "; ")
  }, "")
  out
}

#' Check a coded table of oversight reports
#'
#' Flags the cells of a table coded from oversight reports (by
#' [ovrpar()], by hand, or by any other means) that are not of the form
#' their field calls for, so that a parsing slip or a typing error is
#' found before the table is analysed. Nothing is changed.
#'
#' The form of a column comes from the rule of its field in `schema`:
#' a `date` field must hold real dates as `YYYY-MM-DD`; a `num`, `tag` or
#' `before` field whole numbers from 0; a `flag` field `TRUE`, `FALSE`,
#' yes or no (in any case, or y and n). `shape` adds columns that must
#' follow a pattern, such as case numbers: in a pattern `9` stands for a
#' digit, `A` for a capital letter, `a` for a small letter and any other
#' character for itself, with alternatives separated by `"|"`
#' (`"99-AAA-999|99-AA-999"`). Columns of `schema` that the table lacks
#' are passed over. A missing value (`NA` or `""`) is not an error.
#'
#' @param data The coded table: a data frame, one row per report.
#' @param schema The fields, as [ovrsch()] gives them; only `name` and
#'   `rule` are read.
#' @param shape `NULL`, or patterns named by column: a named character
#'   vector.
#' @return A data frame with one row per row of `data`: `row`, `issues`
#'   (how many cells are flagged) and `what` (each flagged column with its
#'   form, `"date_of_incident_iso: date"`, separated by `"; "`).
#' @examples
#' d <- data.frame(case_number = c("23-OFD-001", "23-OFD-1", NA),
#'                 date_of_incident_iso = c("2023-01-05", "2023-02-30", ""),
#'                 number_of_subject_officials = c(1, 2.5, NA),
#'                 charges_recommended = c(FALSE, NA, TRUE))
#' ovrchk(d, shape = c(case_number = "99-AAA-999|99-AA-999"))
#' @export
ovrchk <- function(data, schema = ovrsch(), shape = NULL) {
  f <- "ovrchk"
  if (!is.data.frame(data)) {
    .ckstop(f, "`data` must be a data frame, not %s", class(data)[1L])
  }
  if (!is.data.frame(schema) || !all(c("name", "rule") %in% names(schema))) {
    .ckstop(f, "`schema` must be a data frame with the columns name and rule")
  }
  form <- c(date = 0L, num = 1L, tag = 1L, before = 1L, flag = 2L)
  nm <- as.character(schema$name)
  k <- form[as.character(schema$rule)]
  use <- !is.na(k) & nm %in% names(data)
  cols <- nm[use]
  kind <- unname(k[use])
  pat <- rep("", length(cols))
  if (!is.null(shape)) {
    .ovrshp(shape, data, f)
    cols <- c(cols, names(shape))
    kind <- c(kind, rep(3L, length(shape)))
    pat <- c(pat, enc2utf8(unname(shape)))
  }
  what <- c("date", "count", "flag", "shape")
  bad <- vapply(seq_along(cols), function(j) {
    .ovrbad(data[[cols[j]]], kind[j], pat[j])
  }, logical(nrow(data)))
  bad <- matrix(bad, nrow(data))
  lab <- paste0(cols, ": ", what[kind + 1L])
  data.frame(
    row = seq_len(nrow(data)),
    issues = as.integer(rowSums(bad)),
    what = vapply(seq_len(nrow(data)), function(i) {
      paste(lab[bad[i, ]], collapse = "; ")
    }, ""),
    stringsAsFactors = FALSE
  )
}

# Text: a character vector (a vector of NA only is text that is missing),
# valid in its encoding, as UTF-8.
.ovrin <- function(x, what, f) {
  if (is.logical(x) && length(x) && all(is.na(x))) {
    x <- as.character(x)
  }
  if (!is.character(x) || !is.null(dim(x))) {
    .ckstop(f, "`%s` must be text (a character vector), not %s", what,
            class(x)[1L])
  }
  enc2utf8(.cktxt(x, what, f))
}

# NA, TRUE or FALSE.
.ovrhtml <- function(html, f) {
  if (!is.logical(html) || length(html) != 1L) {
    .ckstop(f, "`html` must be NA, TRUE or FALSE")
  }
  html
}

# The rules of a schema, in the order of their codes in src/ovr.cpp.
.ovrrules <- c("after", "num", "date", "tag", "sent", "prev", "lines",
               "name", "flag", "before", "word")

# A schema: the columns name, rule, cue, alt and sect, text with no NA,
# known rules and names given once.
.ovrschk <- function(schema, f) {
  cols <- c("name", "rule", "cue", "alt", "sect")
  if (!is.data.frame(schema) || !all(cols %in% names(schema))) {
    .ckstop(f, "`schema` must be a data frame with the columns %s",
            toString(cols))
  }
  s <- lapply(schema[cols], function(v) {
    if (is.factor(v)) as.character(v) else v
  })
  if (!all(vapply(s, is.character, NA)) || anyNA(unlist(s)) ||
        !nrow(schema)) {
    .ckstop(f, "the columns of `schema` must be text with no missing values")
  }
  bad <- setdiff(s$rule, .ovrrules)
  if (length(bad)) {
    .ckstop(f, "`schema` has a rule that is not one of %s: \"%s\"",
            toString(.ovrrules), bad[1L])
  }
  if (!all(nzchar(s$name)) || anyDuplicated(s$name)) {
    .ckstop(f, "the names of `schema` must be given once each")
  }
  s <- lapply(s, enc2utf8)
  s
}

# Words to look for: text with no missing or empty value and no "|", one
# of them when one; joined by "|" for the kernels.
.ovrwords <- function(x, what, f, one = FALSE) {
  .cktxt(x, what, f)
  if (!is.character(x) || !length(x) || anyNA(x) || !all(nzchar(x)) ||
        any(grepl("|", x, fixed = TRUE)) || (one && length(x) != 1L)) {
    .ckstop(f, "`%s` must be %s, not empty, with no \"|\"", what,
            if (one) "one text" else "text")
  }
  paste(enc2utf8(x), collapse = "|")
}

# The keyword groups: a preset, a named list of words, or words.
.ovrvoc <- function(vocab, f) {
  if (is.character(vocab) && length(vocab) == 1L && identical(vocab, "mhr")) {
    return(.Call(C_clab_r_ovrvoc, vocab))
  }
  if (is.list(vocab) && !is.data.frame(vocab)) {
    nm <- names(vocab)
    ok <- length(vocab) && !is.null(nm) && !anyNA(nm) && all(nzchar(nm)) &&
      !anyDuplicated(nm) && all(vapply(vocab, is.character, NA))
    w <- unlist(vocab, use.names = FALSE)
    if (!ok || !length(w) || anyNA(w) || !all(nzchar(w))) {
      .ckstop(f, paste("`vocab` must be \"mhr\", a list of words named by",
                       "group, or words"))
    }
    return(list(group = rep(nm, lengths(vocab)), term = enc2utf8(w)))
  }
  # a word is also a group's name, and R's names hold 10000 bytes at most
  if (!is.character(vocab) || !length(vocab) || anyNA(vocab) ||
        !all(nzchar(vocab)) || any(nchar(vocab, "bytes") > 10000L)) {
    .ckstop(f, paste("`vocab` must be \"mhr\", a list of words named by",
                     "group, or words"))
  }
  list(group = vocab, term = enc2utf8(vocab))
}

# Patterns named by columns of data: text, none missing or empty.
.ovrshp <- function(shape, data, f) {
  if (!is.character(shape) || is.null(names(shape)) || anyNA(shape) ||
        !all(nzchar(shape)) || !all(names(shape) %in% names(data))) {
    .ckstop(f, "`shape` must be patterns named by columns of `data`")
  }
}

# The cells of a column not of their form (kind, pat); a missing value
# is not flagged.
.ovrbad <- function(v, kind, pat) {
  if (is.numeric(v)) {
    # whole numbers in full (as.character() writes 1e+05)
    w <- as.character(v)
    i <- !is.na(v) & v == round(v)
    w[i] <- sprintf("%.0f", v[i])
    w[is.na(v)] <- NA_character_
    v <- w
  }
  v <- enc2utf8(trimws(as.character(v)))
  v[!is.na(v) & !nzchar(v)] <- NA_character_
  ok <- .Call(C_clab_r_ovrchk, v, kind, pat)
  !is.na(ok) & !ok
}
