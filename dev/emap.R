# Where every R expectation of a base revision went: to a C++ check of the
# self-test (src/t_*.cpp, by id) or still in R (file:line of the working
# tree). Run from the package root:
#   Rscript dev/emap.R [base] [out.tsv]       base defaults to main
# One row per expect_*() call of the base's tests/testthat/*.R. Exits 1 if
# any expectation is found in neither place.
#
# Kinds: a value by hand, b enumeration or simulation, c kernel guard by
# .Call, d R-level message, e R-level shape, f rmorie parity, g plot,
# h proof table. They are read off the code (see kind() below); the split
# between a and b, and between a and e for what stays in R, is a heuristic.
# A C++ check's id is the base file's stem and line: disp3, frns2_70,
# crim9 (test-crimlab.R), edge52.
args <- commandArgs(TRUE)
base <- if (length(args) >= 1L) args[[1L]] else "main"
out <- if (length(args) >= 2L) args[[2L]] else ""

norm <- function(x) gsub("\\s+", " ", paste(x, collapse = " "))
txt <- function(e) norm(deparse(e, width.cutoff = 500L))

# The expect_*() calls of one parsed file: line, test name, text, call.
calls <- function(src, file) {
  ex <- parse(text = src, keep.source = TRUE)
  pd <- getParseData(ex)
  ep <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & grepl("^expect_", pd$text), ]
  tt <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "test_that", ]
  tn <- vapply(seq_len(nrow(tt)), function(i) {
    s <- pd[pd$token == "STR_CONST" & pd$line1 >= tt$line1[i], ]
    s <- s[order(s$line1, s$col1), ]
    gsub("^\"|\"$", "", s$text[1L])
  }, "")
  # line ranges of loops: for, and the apply family around an expectation
  lp <- pd[pd$token == "FOR" | (pd$token == "SYMBOL_FUNCTION_CALL" &
    pd$text %in% c("replicate", "apply", "vapply", "sapply", "lapply")), ]
  rng <- t(vapply(lp$parent, function(p) {
    q <- pd[pd$id == p, ]
    if (q$token != "expr") q <- pd[pd$id == q$parent, ]
    c(q$line1, q$line2)
  }, numeric(2)))
  rows <- lapply(seq_len(nrow(ep)), function(i) {
    l <- ep$line1[i]
    e <- pd[pd$id == pd$parent[pd$id == ep$parent[i]], ]
    k <- max(c(0L, which(tt$line1 <= l)))
    data.frame(
      file = file, line = l, test = if (k == 0L) "" else tn[k],
      fn = ep$text[i], text = txt(str2lang(getParseText(pd, e$id))),
      loop = length(rng) > 0L && any(rng[, 1L] <= l & rng[, 2L] >= l),
      skip = FALSE
    )
  })
  r <- do.call(rbind, rows)
  # tests that skip without rmorie
  for (i in seq_len(nrow(tt))) {
    q <- pd[pd$id == pd$parent[pd$id == tt$parent[i]], ]
    body <- getParseText(pd, q$id)
    if (grepl("skip_if_not_installed", body, fixed = TRUE)) {
      r$skip[r$line >= q$line1 & r$line <= q$line2] <- TRUE
    }
  }
  r
}

# The table-driven error cases of the working tree, rebuilt as the call
# they stand for: errs("f", list(list(msg, args...), ...)).
errcases <- function(path) {
  pd <- getParseData(parse(path, keep.source = TRUE))
  head <- function(id) {  # the function a call expression calls, or ""
    k <- pd[pd$parent == id & pd$token == "expr", ]
    k <- k[order(k$line1, k$col1), ][1L, ]
    f <- pd[pd$parent == k$id & pd$token == "SYMBOL_FUNCTION_CALL", ]
    if (nrow(f)) f$text else ""
  }
  up <- function(id) pd$parent[pd$id == id]
  li <- pd[pd$token == "SYMBOL_FUNCTION_CALL" & pd$text == "list", ]
  res <- NULL
  for (j in seq_len(nrow(li))) {
    case <- up(up(li$id[j]))
    tab <- up(case)
    if (head(tab) != "list" || head(up(tab)) != "errs") next
    e <- str2lang(getParseText(pd, case))
    fe <- str2lang(getParseText(pd, up(tab)))[[2L]]
    a <- as.list(e)[-1L]
    cl <- as.call(c(as.name(fe), a[-1L]))
    res <- rbind(res, data.frame(text = txt(call("expect_error", cl, a[[1L]])),
                                 line = li$line1[j]))
  }
  res
}

kind <- function(r, moved) {
  t <- r$text
  if (r$skip) return("f")
  if (grepl("proof", r$test) || r$file == "test-proof.R") return("h")
  if (grepl("tipplt|eceplt", t) || grepl("plot", r$test)) return("g")
  if (grepl(".Call(", t, fixed = TRUE)) return("c")
  if (r$fn == "expect_error") return("d")
  if (moved) {
    b <- r$loop || grepl("brute|enumerat|extremes|grid|any data|monoton|add up",
                         r$test)
    return(if (b) "b" else "a")
  }
  shape <- r$fn %in% c("expect_identical", "expect_named", "expect_false",
                       "expect_no_warning") ||
    grepl("names\\(|is\\.na\\(|\\[\\[\"(group|term|locus|off|vic|part)\"\\]\\]",
          t) || r$file == "test-chk.R"
  if (shape) "e" else "a"
}

stem <- function(f) {
  s <- sub("^test-(.*)\\.R$", "\\1", f)
  if (s == "crimlab") "crim" else if (grepl("[0-9]$", s)) paste0(s, "_") else s
}

# the base's expectations
fs <- system2("git", c("ls-tree", "--name-only", base, "tests/testthat/"),
              stdout = TRUE)
fs <- fs[grepl("/test-.*\\.R$", fs)]
orig <- do.call(rbind, lapply(fs, function(f) {
  calls(system2("git", c("show", paste0(base, ":", f)), stdout = TRUE),
        basename(f))
}))

# ids of the C++ checks
cpp <- unlist(lapply(Sys.glob("src/t_*.cpp"), readLines))
ids <- unique(unlist(regmatches(cpp, gregexpr("(?<=, \")[a-z0-9_]+(?=\")", cpp,
                                              perl = TRUE))))

# the working tree's R expectations, plain and table-driven
now <- NULL
for (f in Sys.glob("tests/testthat/test-*.R")) {
  c1 <- calls(readLines(f), basename(f))
  now <- rbind(now, data.frame(file = basename(f), line = c1$line,
                               text = c1$text))
  e1 <- errcases(f)
  if (!is.null(e1)) {
    now <- rbind(now, data.frame(file = basename(f), line = e1$line,
                                 text = e1$text))
  }
}
now$used <- FALSE

orig$id <- ""
orig$kept <- ""
orig$kind <- ""
for (i in seq_len(nrow(orig))) {
  id <- paste0(stem(orig$file[i]), orig$line[i])
  if (id %in% ids) orig$id[i] <- id
  m <- which(now$text == orig$text[i] & !now$used)
  # the same file first
  m <- m[order(now$file[m] != orig$file[i])]
  if (length(m)) {
    now$used[m[1L]] <- TRUE
    orig$kept[i] <- paste0(now$file[m[1L]], ":", now$line[m[1L]])
  }
  orig$kind[i] <- kind(orig[i, ], nzchar(orig$id[i]) && !nzchar(orig$kept[i]))
}
orig$to <- ifelse(nzchar(orig$id) & nzchar(orig$kept),
                  paste0("c++ ", orig$id, "; kept in R: ", orig$kept),
                  ifelse(nzchar(orig$id), paste0("c++ ", orig$id),
                         ifelse(nzchar(orig$kept),
                                paste0("kept in R: ", orig$kept), "MISSING")))
tab <- orig[, c("file", "line", "test", "kind", "to")]
write.table(tab, if (nzchar(out)) out else stdout(), sep = "\t",
            quote = FALSE, row.names = FALSE)
moved <- nzchar(orig$id) & !nzchar(orig$kept)
message(sprintf("%d expectations: %d to C++, %d kept in R (%d also in C++), %d missing",
                nrow(orig), sum(moved), sum(nzchar(orig$kept)),
                sum(nzchar(orig$id) & nzchar(orig$kept)),
                sum(orig$to == "MISSING")))
message("by kind: ", paste(names(table(orig$kind)), table(orig$kind),
                           sep = "=", collapse = " "))
message("moved by kind: ", paste(names(table(orig$kind[moved])),
                                 table(orig$kind[moved]), sep = "=",
                                 collapse = " "))
quit(status = if (any(orig$to == "MISSING")) 1L else 0L)
