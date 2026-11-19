# A change table rendered to a self-contained HTML page or to a PDF
# written by the package itself (no device, no external program), and
# the writer that picks the format from the file name. Colour encodes the
# DIRECTION of an improvement, not the sign of a change: with `better`
# set, a rise in days in segregation or in reported assaults is drawn as
# worse and a rise in a completion rate as better.

# the palettes are .YOY_PALETTES (R/plotx.R), shared with the plots

#' Render a change table to HTML or PDF
#'
#' `yoyhtml()` writes one self-contained HTML page (every style inlined,
#' nothing fetched, so the file renders the same years later): the table,
#' a card per series with its change over the span ([yoysmry()]) and the
#' notes. `yoypdf()` writes a PDF of the table, paginated when it does not
#' fit, with the standard Helvetica font and no other program: crimlab
#' writes the file itself.
#'
#' Colour encodes the direction of an improvement, not the sign of the
#' change: with `better = "lower"` (fewer assaults, fewer days in
#' segregation) a rise is drawn in the colour for worse, with `"higher"`
#' (a completion rate) in the colour for better, and with `"neutral"` (the
#' default) every change in the neutral colour. The palette `"safe"`
#' (blue and orange) stays apart in the common forms of colour blindness;
#' `"mono"` is for print. The first period of each series, with nothing to
#' compare with, is said once under the table rather than on every first
#' row. The PDF holds ASCII only: dashes and arrows become `-` and `+`,
#' other characters their transliteration ([toascii()]).
#'
#' @inheritParams yoylbl
#' @param file The file path. A name with no suffix is given the
#'   format's own (`.html`, `.pdf`).
#' @param title The heading (and the page title).
#' @param sub The line under the heading.
#' @param pal The palette: `"diverging"`, `"safe"` or `"mono"`.
#' @param better Which direction is better: `"neutral"`, `"higher"` or
#'   `"lower"`.
#' @param digits The decimals of the change shown.
#' @param bars `TRUE` draws a bar under each change, proportional to its
#'   size.
#' @param intv `TRUE` shows the interval column of a count table.
#' @param notes `NULL`, or footnotes.
#' @param lag The periods per step of the table, for the summary and the
#'   note on the first periods.
#' @param width,height The PDF page in inches, each from 2 to 200.
#' @return The path written, invisibly.
#' @seealso [yoywr()], [yoycsv()].
#' @examples
#' # simulated police-reported counts by offence and year, any agency
#' d <- data.frame(year = rep(2021:2023, 2),
#'                 offence = rep(c("Assault", "Robbery"), each = 3),
#'                 n = c(812, 845, 790, 120, 131, 98))
#' y <- yoypx(d$n, d$year, d$offence)
#' f <- tempfile(fileext = ".html")
#' yoyhtml(y, f, title = "Reported offences (simulated)", better = "lower")
#' g <- tempfile(fileext = ".pdf")
#' yoypdf(y, g, title = "Reported offences (simulated)")
#' # simulated days in segregation by fiscal year
#' s <- yoylbl(yoypx(c(4100, 3900, 4400), 2021:2023), fylbl)
#' yoyhtml(s, f, title = "Days in segregation (simulated)", better = "lower")
#' unlink(c(f, g))
#' @name yoyhtml
#' @export
yoyhtml <- function(x, file, title = "Year-over-year change", sub = "",
                    pal = "diverging", better = "neutral", digits = 1L,
                    bars = TRUE, intv = TRUE, notes = NULL, lag = 1) {
  d <- .ytab(x, "yoyhtml")
  file <- .cksfx(.y1file(file, "yoyhtml"), "html")
  .yargs(title, sub, notes, "yoyhtml")
  p <- .YOY_PALETTES[[.ypick(pal, names(.YOY_PALETTES), "pal", "yoyhtml")]]
  better <- .ypick(better, c("neutral", "higher", "lower"), "better",
                   "yoyhtml")
  digits <- .ydig(digits, "yoyhtml")
  bars <- .yflag(bars, "bars", "yoyhtml")
  intv <- .yflag(intv, "intv", "yoyhtml")
  s <- .yshow(d, digits, intv, md = FALSE)
  chg <- d[[.ychg(d)]]
  cls <- .ycls(chg, better)
  sc <- suppressWarnings(max(abs(chg), na.rm = TRUE))
  if (!is.finite(sc) || sc == 0) sc <- 1
  rows <- vapply(seq_len(nrow(d)), .yhrow, "", s = s, chg = chg, cls = cls,
                 sc = sc, bars = bars)
  sm <- .clas("yoyhtml", yoysmry(d, lag))
  cards <- vapply(seq_len(nrow(sm)), .yhcard, "", sm = sm)
  html <- .yhpage(title, sub, p, names(s$cols), rows, cards,
                  c(.yfirst(d, lag), notes))
  writeBin(charToRaw(enc2utf8(paste0(paste(html, collapse = "\n"), "\n"))),
           file)
  invisible(file)
}

# Helvetica widths (1/1000 em) of the printable ASCII characters
.yhw <- c(278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278,
          333, 278, 278, rep(556, 10), 278, 278, 584, 584, 584, 556, 1015,
          667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833,
          722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611,
          278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556,
          556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278,
          556, 500, 722, 500, 500, 500, 334, 260, 334, 584)

#' @rdname yoyhtml
#' @export
yoypdf <- function(x, file, title = "Year-over-year change", sub = "",
                   pal = "diverging", better = "neutral", digits = 1L,
                   intv = TRUE, notes = NULL, width = 11, height = 8.5,
                   lag = 1) {
  d <- .ytab(x, "yoypdf")
  file <- .cksfx(.y1file(file, "yoypdf"), "pdf")
  .yargs(title, sub, notes, "yoypdf")
  .ylag(lag, "yoypdf")
  p <- .YOY_PALETTES[[.ypick(pal, names(.YOY_PALETTES), "pal", "yoypdf")]]
  better <- .ypick(better, c("neutral", "higher", "lower"), "better",
                   "yoypdf")
  digits <- .ydig(digits, "yoypdf")
  intv <- .yflag(intv, "intv", "yoypdf")
  .ypsize(width, height)
  s <- .yshow(d, digits, intv, md = FALSE)
  lay <- .yplay(d, s, width, height)
  lay$body <- lapply(s$cols, .yasc)
  lay$cls <- .ycls(d[[.ychg(d)]], better)
  lay$nt <- .yasc(c(.yfirst(d, lag), notes))
  lay$p <- p
  lay$dim <- .yrgb(p$flat)
  lay$title <- .yasc(title)
  lay$sub <- .yasc(sub)
  streams <- vapply(seq_along(lay$pages), .ypstr, "", lay = lay)
  .ypdfw(streams, lay$W, lay$H, .yasc(title), file)
}

#' Write a change table in the format its file name says
#'
#' Picks the writer from the extension of `file` (`.html` or `.htm`,
#' `.pdf`, `.csv`, `.tsv` or `.tab`, `.json`, `.md` or `.markdown`), or
#' from `format`, and passes the other arguments to it.
#'
#' @inheritParams yoylbl
#' @inheritParams yoyhtml
#' @param format `"auto"` (from the extension), or one of `"html"`, `"pdf"`,
#'   `"csv"`, `"tsv"`, `"json"`, `"markdown"`.
#' @param ... The arguments of the chosen writer.
#' @return The path written, invisibly.
#' @seealso [yoycsv()], [yoyhtml()].
#' @examples
#' # simulated police-reported counts by offence and year
#' y <- yoypx(c(812, 845, 790), 2021:2023)
#' f <- tempfile(fileext = ".csv")
#' yoywr(y, f)
#' readLines(f)
#' unlink(f)
#' @export
yoywr <- function(x, file, format = "auto", ...) {
  .y1file(file, "yoywr")
  format <- tolower(.ytxt1(format, "format", "yoywr"))
  if (format == "auto") {
    map <- c(html = "html", htm = "html", pdf = "pdf", csv = "csv",
             tsv = "tsv", tab = "tsv", json = "json", md = "markdown",
             markdown = "markdown")
    ext <- tolower(sub(".*\\.", "", basename(file)))
    if (!ext %in% names(map)) {
      stop("yoywr: cannot tell the format from \"", basename(file),
           "\"; pass `format`.", call. = FALSE)
    }
    format <- map[[ext]]
  }
  f <- switch(format, html = yoyhtml, pdf = yoypdf, csv = yoycsv,
              tsv = yoytsv, json = yoyjs, markdown = yoymd)
  if (is.null(f)) {
    stop("yoywr: unknown format: ", format, ".", call. = FALSE)
  }
  invisible(.clas("yoywr", f(x, file, ...)))
}

#' `x` checked as one of the options `opt`.
#' @noRd
.ypick <- function(x, opt, nm, fn) {
  if (!is.character(x) || length(x) != 1L || !x %in% opt) {
    stop(fn, ": `", nm, "` must be one of ", toString(dQuote(opt, FALSE)),
         ".", call. = FALSE)
  }
  x
}

#' the class of each row: bad, good or flat
#' @noRd
.ycls <- function(chg, better) {
  out <- rep("flat", length(chg))
  if (better == "neutral") {
    return(out)
  }
  up <- !is.na(chg) & chg > 0
  dn <- !is.na(chg) & chg < 0
  out[up] <- if (better == "higher") "good" else "bad"
  out[dn] <- if (better == "higher") "bad" else "good"
  out
}

#' Text escaped for HTML.
#' @noRd
.yesc <- function(s) {
  s <- as.character(s)
  s[is.na(s)] <- ""
  s <- .yctl(s, ws = TRUE)
  s <- gsub("&", "&amp;", s, fixed = TRUE)
  s <- gsub("<", "&lt;", s, fixed = TRUE)
  s <- gsub(">", "&gt;", s, fixed = TRUE)
  s <- gsub("\"", "&quot;", s, fixed = TRUE)
  gsub("'", "&#39;", s, fixed = TRUE)
}

#' The note that the first period has nothing to compare with, NULL when no
#'   row is flagged so.
#' @noRd
.yfirst <- function(d, lag) {
  if (!("flag" %in% names(d)) ||
        !any(!is.na(d$flag) & d$flag == "no earlier period")) {
    return(NULL)
  }
  sprintf(
    "Note: the first %s%s %s nothing to compare with, so %s no change.",
    if (lag == 1L) "period" else sprintf("%s periods", format(lag)),
    if ("group" %in% names(d)) " of each group" else "",
    if (lag == 1L) "has" else "have",
    if (lag == 1L) "it shows" else "they show"
  )
}

#' The style sheet of the HTML table for the palette `p`.
#' @noRd
.ycss <- function(p) paste0(.ycssp(p), .ycsst())

#' the palette's colours and the page's, light and dark
#' @noRd
.ycssp <- function(p) {
  paste0(
    ":root{--bad:", p$bad, ";--good:", p$good, ";--flat:", p$flat,
    ";--badbg:", p$badbg, ";--goodbg:", p$goodbg,
    ";--ink:#1a1a1a;--dim:#5c5c5c;--rule:#e0ddd8;--ground:#faf9f7;",
    "--panel:#ffffff}",
    "@media (prefers-color-scheme:dark){:root:not([data-theme=light]){",
    "--ink:#eceae6;--dim:#a6a29b;--rule:#332f2b;--ground:#171513;",
    "--panel:#201d1a;--badbg:#331f1a;--goodbg:#142a26}}"
  )
}

#' the page's layout and the table's
#' @noRd
.ycsst <- function() {
  paste0(
    "*{box-sizing:border-box}",
    "body{margin:0;background:var(--ground);color:var(--ink);",
    "font:14px/1.5 ui-sans-serif,-apple-system,'Segoe UI',Roboto,",
    "'Helvetica Neue',Arial,sans-serif}",
    "main{max-width:1000px;margin:0 auto;padding:2.5rem 1.25rem 4rem}",
    "h1{font-size:1.6rem;margin:0 0 .25rem;letter-spacing:-.01em}",
    "h2{font-size:1rem;text-transform:uppercase;letter-spacing:.08em;",
    "color:var(--dim);margin:2.5rem 0 .75rem;font-weight:600}",
    "p.sub{color:var(--dim);margin:0 0 1.75rem}",
    ".scroll{overflow-x:auto;border:1px solid var(--rule);border-radius:6px;",
    "background:var(--panel)}",
    "table{border-collapse:collapse;width:100%;min-width:520px}",
    "th{text-align:left;font-size:.72rem;text-transform:uppercase;",
    "letter-spacing:.06em;color:var(--dim);font-weight:600;",
    "padding:.65rem .75rem;border-bottom:1px solid var(--rule);",
    "white-space:nowrap}",
    "td{padding:.55rem .75rem;border-bottom:1px solid var(--rule);",
    "vertical-align:middle}",
    "tbody tr:last-child td{border-bottom:0}",
    "td.num,td.ci{text-align:right;font-variant-numeric:tabular-nums}",
    "td.per{font-variant-numeric:tabular-nums;white-space:nowrap}",
    "td.grp{font-weight:600}",
    "td.ci{color:var(--dim);font-size:.85rem;white-space:nowrap}",
    "td.pct{position:relative;text-align:right;white-space:nowrap;",
    "font-variant-numeric:tabular-nums}",
    "td.pct .arw{margin-right:.3rem;font-size:.8rem}",
    "td.pct.bad{color:var(--bad)}td.pct.good{color:var(--good)}",
    "td.pct.flat{color:var(--flat)}",
    ".bar{position:absolute;left:0;bottom:0;height:2px;opacity:.55}",
    ".bar.bad{background:var(--bad)}.bar.good{background:var(--good)}",
    ".bar.flat{background:var(--flat)}",
    "tr.bad td.pct{background:var(--badbg)}",
    "tr.good td.pct{background:var(--goodbg)}",
    "tr.note td{border-bottom:1px solid var(--rule);color:var(--dim);",
    "font-size:.8rem;padding:.2rem .75rem .5rem}",
    ".cards{display:grid;gap:.75rem;",
    "grid-template-columns:repeat(auto-fit,minmax(180px,1fr))}",
    ".card{border:1px solid var(--rule);border-radius:6px;padding:.9rem;",
    "background:var(--panel)}",
    ".card .lab{font-size:.72rem;text-transform:uppercase;",
    "letter-spacing:.06em;color:var(--dim);font-weight:600}",
    ".card .big{font-size:1.5rem;margin:.2rem 0 .3rem;",
    "font-variant-numeric:tabular-nums}",
    ".card .cap{color:var(--dim);font-size:.8rem}",
    "ul.notes{color:var(--dim);font-size:.85rem;padding-left:1.1rem}",
    "p.gen{color:var(--dim);font-size:.75rem;margin-top:2.5rem}"
  )
}

#' `title`, `sub` and `notes` checked.
#' @noRd
.yargs <- function(title, sub, notes, fn) {
  .ytxt1(title, "title", fn)
  .ytxt1(sub, "sub", fn)
  if (!is.null(notes) && (!is.character(notes) || anyNA(notes))) {
    stop(fn, ": `notes` must be text (a character vector), or NULL.",
         call. = FALSE)
  }
  .cktxt(notes, "notes", fn)
}

#' the PDF's text: dashes and arrows as ASCII, then toascii()
#' @noRd
.yasc <- function(x) {
  out <- as.character(x)
  for (p in list(c("\u2014", "-"), c("\u2013", "-"), c("\u25b2", "+"),
                 c("\u25bc", "-"), c("\u00b7", "|"), c("\u2500", "-"),
                 c("\u21b3", ">"))) {
    out <- gsub(p[1L], p[2L], out, fixed = TRUE)
  }
  toascii(out)
}

#' The width of each string at a font size, from the table of character
#'   widths.
#' @noRd
.ywid <- function(s, size) {
  vapply(s, function(t) {
    k <- utf8ToInt(t) - 31L
    sum(.yhw[k[k >= 1L & k <= 95L]]) * size / 1000
  }, 0, USE.NAMES = FALSE)
}

#' Text escaped for a PDF string.
#' @noRd
.ypdfs <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("(", "\\(", s, fixed = TRUE)
  gsub(")", "\\)", s, fixed = TRUE)
}

#' A hex colour as three numbers from 0 to 1, for PDF.
#' @noRd
.yrgb <- function(hex) {
  v <- strtoi(substring(hex, c(2L, 4L, 6L), c(3L, 5L, 7L)), 16L) / 255
  paste(formatC(v, format = "f", digits = 3), collapse = " ")
}

#' yoypdf(): the page, each side a number from 2 to 200 inches
#' @noRd
.ypsize <- function(width, height) {
  for (a in list(list("width", width), list("height", height))) {
    v <- a[[2L]]
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v < 2 ||
          v > 200) {
      stop("yoypdf: `", a[[1L]], "` must be one number from 2 to 200.",
           call. = FALSE)
    }
  }
}

#' The layout of yoypdf(): page size and margin in points, the edges of
#' the columns, the row height and the rows of each page
#' @noRd
.yplay <- function(d, s, width, height) {
  hdr <- .yasc(names(s$cols))
  W <- 72 * width
  H <- 72 * height
  m <- 0.6 * 72 / 6
  nc <- length(hdr)
  numfrom <- as.integer("group" %in% names(d)) + 2L
  w <- rep(1, nc)
  w[seq_len(numfrom - 1L)] <- 1.3
  per <- max(5L, as.integer(floor((0.82 * height - 0.22 * 1.1) / 0.22)))
  pages <- split(seq_len(nrow(d)), ceiling(seq_len(nrow(d)) / per))
  if (!length(pages)) pages <- list(integer())
  list(hdr = hdr, W = W, H = H, m = m, nc = nc, numfrom = numfrom,
       xs = m + (W - 2 * m) * cumsum(c(0, w / sum(w))), rowh = 0.22 * 72,
       pages = pages)
}

#' yoypdf(): a coordinate to the hundredth of a point, a half going up. A
#' right-aligned text ends on thousandths of a point, so a half is common,
#' and "%.2f" alone sends it up or down by the last bit of the sum, which
#' differs between machines: the file must be the same bytes on each.
#' @noRd
.ypc <- function(x) floor(x * 100 + 0.5 + 1e-6) / 100

#' yoypdf(): a line of text at (x0, y0), its right end there when `right`
#' @noRd
.yptxt <- function(x0, y0, s, size, font = "F1", rgb = "0 0 0",
                   right = FALSE) {
  if (right) x0 <- x0 - .ywid(s, size)
  sprintf("BT /%s %g Tf %s rg %.2f %.2f Td (%s) Tj ET", font, size, rgb,
          .ypc(x0), .ypc(y0), .ypdfs(s))
}

#' yoypdf(): the content stream of page k: heading, column heads, rows,
#' page number, and the notes on the last page
#' @noRd
.ypstr <- function(k, lay) {
  m <- lay$m
  W <- lay$W
  H <- lay$H
  dim <- lay$dim
  np <- length(lay$pages)
  top <- H - m - 0.12 * H
  rule <- function(y) {
    sprintf("%s RG 0.7 w %.2f %.2f m %.2f %.2f l S", dim, .ypc(m), .ypc(y),
            .ypc(W - m), .ypc(y))
  }
  out <- c(.yptxt(m, H - m - 18, lay$title, 18, "F2"),
           .yptxt(m, H - m - 32, lay$sub, 9, rgb = dim),
           rule(top + 10))
  for (j in seq_len(lay$nc)) {
    r <- j >= lay$numfrom && lay$hdr[j] != "note"
    out <- c(out, .yptxt(if (r) lay$xs[j + 1L] - 3 else lay$xs[j] + 3, top,
                         lay$hdr[j], 8, "F2", dim, r))
  }
  out <- c(out, rule(top - 5), .yprows(lay, lay$pages[[k]], top))
  if (np > 1L) {
    out <- c(out, .yptxt(W - m, m, sprintf("page %d of %d", k, np), 7,
                         rgb = dim, right = TRUE))
  }
  if (k == np && length(lay$nt)) {
    yy <- m + 10 + 10 * length(lay$nt)
    for (q in lay$nt) {
      out <- c(out, .yptxt(m, yy, q, 7.5, rgb = dim))
      yy <- yy - 10
    }
  }
  paste(out, collapse = "\n")
}

#' yoypdf(): the table rows `rows` of a page, from under the column heads
#' @noRd
.yprows <- function(lay, rows, top) {
  out <- character(0)
  y <- top - lay$rowh * 1.1
  pc <- which(lay$hdr %in% c("change %", "points"))
  for (i in rows) {
    for (j in seq_len(lay$nc)) {
      r <- j >= lay$numfrom
      col <- if (j == pc) .yrgb(lay$p[[lay$cls[i]]]) else "0 0 0"
      out <- c(out, .yptxt(if (r) lay$xs[j + 1L] - 3 else lay$xs[j] + 3, y,
                           lay$body[[j]][i], 9, rgb = col, right = r))
    }
    y <- y - lay$rowh
  }
  out
}

#' the PDF file: catalog, pages, two standard fonts, a page and a content
#' stream per page, the cross-reference table
#' @noRd
.ypdfw <- function(streams, W, H, title, file) {
  np <- length(streams)
  kids <- 5L + 2L * (seq_len(np) - 1L)
  obj <- c(
    "<< /Type /Catalog /Pages 2 0 R >>",
    sprintf("<< /Type /Pages /Kids [%s] /Count %d >>",
            paste(kids, "0 R", collapse = " "), np),
    paste("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica",
          "/Encoding /WinAnsiEncoding >>"),
    paste("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica-Bold",
          "/Encoding /WinAnsiEncoding >>")
  )
  for (k in seq_len(np)) {
    obj <- c(obj, sprintf(paste0(
      "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] /Resources ",
      "<< /Font << /F1 3 0 R /F2 4 0 R >> >> /Contents %d 0 R >>"
    ), W, H, kids[k] + 1L),
    sprintf("<< /Length %d >>\nstream\n%s\nendstream",
            nchar(streams[k], type = "bytes") + 0L, streams[k]))
  }
  obj <- c(obj, sprintf("<< /Title (%s) /Producer (crimlab) >>",
                        .ypdfs(title)))
  head <- "%PDF-1.4\n"
  parts <- sprintf("%d 0 obj\n%s\nendobj\n", seq_along(obj), obj)
  off <- nchar(head, type = "bytes") +
    cumsum(c(0, nchar(parts, type = "bytes")))[seq_along(parts)]
  xref <- paste0("xref\n0 ", length(obj) + 1L, "\n0000000000 65535 f \n",
                 paste(sprintf("%010d 00000 n \n", off), collapse = ""))
  at <- nchar(head, type = "bytes") + sum(nchar(parts, type = "bytes"))
  tr <- sprintf(paste0("trailer\n<< /Size %d /Root 1 0 R /Info %d 0 R >>\n",
                       "startxref\n%d\n%%%%EOF\n"),
                length(obj) + 1L, length(obj), at)
  writeBin(charToRaw(paste0(head, paste(parts, collapse = ""), xref, tr)),
           file)
  invisible(file)
}

#' yoyhtml(), yoypdf(), yoywr(): `file` is one path
#' @noRd
.y1file <- function(file, fn) {
  if (!is.character(file) || length(file) != 1L || is.na(file)) {
    stop(fn, ": `file` must be one file path.", call. = FALSE)
  }
  .ckout(file, fn)
}

#' yoyhtml(): row r of the table, with its note row when it has a note
#' @noRd
.yhrow <- function(r, s, chg, cls, sc, bars) {
  hc <- names(s$cols)
  pc <- which(hc %in% c("change %", "points"))
  cells <- vapply(seq_along(hc), function(j) {
    v <- s$cols[[j]][r]
    if (j == pc) {
      return(.yhpct(v, chg[r], cls[r], sc, bars))
    }
    k <- switch(hc[j], group = "grp", prd = "per", interval = "ci", "num")
    paste0("<td class=\"", k, "\">", if (k %in% c("grp", "per")) .yesc(v)
           else v, "</td>")
  }, "")
  f <- s$flag[r]
  note <- if (is.na(f) || f == "no earlier period") "" else
    sprintf("<tr class=\"note\"><td colspan=\"%d\">note: %s</td></tr>",
            length(hc), .yesc(f))
  paste0("<tr class=\"", cls[r], "\">", paste(cells, collapse = ""),
         "</tr>", note)
}

#' yoyhtml(): the cell of a change: its arrow, its value and its bar
#' @noRd
.yhpct <- function(v, ch, cl, sc, bars) {
  a <- if (is.na(ch)) "" else if (ch > 0) "\u25b2" else
    if (ch < 0) "\u25bc" else "\u2013"
  bar <- if (bars && !is.na(ch) && is.finite(ch)) {
    sprintf("<span class=\"bar %s\" style=\"width:%.1f%%\"></span>",
            cl, min(100, 100 * abs(ch) / sc))
  } else {
    ""
  }
  paste0("<td class=\"pct ", cl, "\"><span class=\"arw\">",
         a, "</span><span class=\"val\">", v, "</span>", bar,
         "</td>")
}

#' yoyhtml(): the card of row r of the summary
#' @noRd
.yhcard <- function(r, sm) {
  paste0(
    "<div class=\"card\"><div class=\"lab\">",
    .yesc(if ("group" %in% names(sm)) sm$group[r] else "overall"),
    "</div><div class=\"big\">", .ypct(sm$totalpct[r], 1L, FALSE),
    "</div><div class=\"cap\">",
    .yesc(sprintf("%s to %s", sm$from[r], sm$to[r])),
    "</div><div class=\"cap\">",
    .yesc(sprintf("compound %s per period",
                  .ypct(sm$cagrpct[r], 1L, FALSE))),
    "</div><div class=\"cap\">",
    .yesc(sprintf("%d up, %d down, %d unchanged", sm$up[r], sm$down[r],
                  sm$flat[r])),
    "</div></div>")
}

#' yoyhtml(): the page's lines
#' @noRd
.yhpage <- function(title, sub, p, hc, rows, cards, nt) {
  c(
    "<!DOCTYPE html>",
    "<html lang=\"en\"><head>",
    "<meta charset=\"utf-8\">",
    "<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">",
    paste0("<title>", .yesc(title), "</title>"),
    paste0("<style>", .ycss(p), "</style>"),
    "</head><body><main>",
    paste0("<h1>", .yesc(title), "</h1>"),
    paste0("<p class=\"sub\">", .yesc(sub), "</p>"),
    "<div class=\"scroll\"><table><thead><tr>",
    paste(paste0("<th>", .yesc(hc), "</th>"), collapse = ""),
    "</tr></thead><tbody>",
    rows,
    "</tbody></table></div>",
    paste0("<h2>Across the span</h2><div class=\"cards\">",
           paste(cards, collapse = ""), "</div>"),
    if (length(nt)) {
      c("<ul class=\"notes\">", paste0("<li>", .yesc(nt), "</li>"), "</ul>")
    },
    paste0("<p class=\"gen\">",
           .yesc(sprintf("crimlab %s", utils::packageVersion("crimlab"))),
           "</p>"),
    "</main></body></html>"
  )
}
