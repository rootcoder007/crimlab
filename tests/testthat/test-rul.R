# rul(), valrules() and the ready-made rules (R/rul.R). The command line
# and the other languages run the same rules in C++ (cli/jval.inc); the
# cases of cli/test.tsv hold the two to the same answers.

d <- data.frame(
  year = c(2019, 2020, 2021, 2021),
  n = c(410, 395, 9000, 402),
  div = c("11", "14", "14", "5X"),
  m = c(1, NA, 3, 4),
  stringsAsFactors = FALSE
)

test_that("rul() declares a rule and valrules() reports the rows", {
  r <- rul(
    "n_small",
    function(v) v < 1000,
    col = "n",
    msg = "a count is too large"
  )
  expect_s3_class(r, "clrule")
  expect_named(r, c("rule", "name", "col", "sev", "msg", "expr"))
  expect_identical(r[["rule"]], "fn")
  expect_identical(r[["sev"]], "warning")
  v <- valrules(d, r)
  expect_s3_class(v, "valrules")
  expect_identical(attr(v, "nrule"), 1L)
  expect_identical(
    v[["n_small"]],
    list(
      severity = "warning",
      message = "a count is too large (1 of 4 rows)",
      nfail = 1L,
      rows = 3L
    )
  )
  # a rule on the table, with the default message and a missing answer
  t <- rul("m_known", function(x) x[["m"]] > 0, sev = "fatal")
  expect_null(t[["col"]])
  expect_identical(t[["msg"]], "Rule 'm_known' failed")
  expect_identical(valrules(d, t)[["m_known"]][["rows"]], 2L)
  # one answer for the whole table
  expect_identical(
    valrules(d, rul("few", function(x) nrow(x) < 3))[["few"]],
    list(severity = "warning", message = "Rule 'few' failed")
  )
  expect_length(valrules(d, rul("some", function(x) nrow(x) > 3)), 0L)
  # more than ten failing rows: the first ten
  big <- data.frame(x = 1:30)
  w <- valrules(big, rul("neg", function(v) v < 0, col = "x"))[["neg"]]
  expect_identical(w[["nfail"]], 30L)
  expect_identical(w[["rows"]], 1:10)
})

test_that("a broken rule is reported, never passed", {
  v <- valrules(
    d,
    list(
      rul("stops", function(v) stop("no such thing"), col = "n"),
      rul("number", function(v) 1, col = "n"),
      rul("empty", function(v) logical(0), col = "n"),
      rul("absent", function(v) FALSE, col = "nocol")
    )
  )
  expect_named(v, c("stops", "number", "empty"))
  expect_match(
    v[["stops"]][["message"]],
    "could not be evaluated: no such thing",
    fixed = TRUE
  )
  expect_match(
    v[["number"]][["message"]],
    "the rule returned numeric, not TRUE or FALSE",
    fixed = TRUE
  )
  expect_match(
    v[["empty"]][["message"]],
    "the rule returned logical, not TRUE or FALSE",
    fixed = TRUE
  )
  # a table of no rows passes a row-wise rule
  expect_length(valrules(d[0, ], rul("r", function(v) v > 0, col = "n")), 0L)
  # on one row, an answer is that row's
  expect_identical(
    valrules(d[3, ], rul("r", function(v) v < 1000, col = "n"))[["r"]][[
      "rows"
    ]],
    1L
  )
})

test_that("the ready-made rules find what they are for", {
  v <- valrules(
    d,
    list(
      rulebtw("year", 2000, 2020),
      ruleincr("year", strict = TRUE),
      rulemad("n", 5),
      rulerx("div", "^[0-9]{2}$"),
      ruledstn(c("year", "div")),
      rulecolcnt(3),
      rulecmpl(),
      rulenull("m"),
      ruleuniq("year"),
      ruleinset("div", c("11", "14")),
      ruledstn("year"),
      ruleincr("year"),
      rulecolcnt(4),
      ruledstn(),
      ruledstn("nocol")
    )
  )
  expect_identical(attr(v, "nrule"), 15L)
  expect_named(
    v,
    c(
      "year_between",
      "year_increasing",
      "n_within_mads",
      "div_regex",
      "col_count",
      "complete_rows",
      "m_not_null",
      "year_unique",
      "div_in_set",
      "distinct_rows"
    )
  )
  expect_identical(v[["year_between"]][["rows"]], 3:4)
  expect_identical(
    v[["year_increasing"]],
    list(
      severity = "warning",
      message = "Column 'year' is not strictly increasing"
    )
  )
  expect_identical(v[["n_within_mads"]][["rows"]], 3L)
  expect_identical(v[["div_regex"]][["rows"]], 4L)
  expect_identical(v[["col_count"]][["message"]], "Expected 3 columns")
  expect_identical(v[["complete_rows"]][["rows"]], 2L)
  expect_identical(v[["m_not_null"]][["rows"]], 2L)
  expect_identical(v[["year_unique"]][["rows"]], 4L)
  expect_identical(v[["div_in_set"]][["rows"]], 4L)
  # two rules of one name: the later one stays (rows repeated on year)
  expect_identical(v[["distinct_rows"]][["rows"]], 4L)
  expect_identical(
    v[["distinct_rows"]][["message"]],
    "Some rows are duplicated on (year) (1 of 4 rows)"
  )
})

test_that("missing values pass a rule about values unless told otherwise", {
  rows <- function(r) {
    v <- valrules(d, r)
    if (length(v)) v[[1L]][["rows"]]
  }
  expect_null(rows(rulebtw("m", 0, 10)))
  expect_identical(rows(rulebtw("m", 0, 10, napass = FALSE)), 2L)
  expect_null(rows(ruleinset("m", c(1, 3, 4))))
  expect_identical(rows(ruleinset("m", c(1, 3, 4), napass = FALSE)), 2L)
  expect_null(rows(rulerx("m", "^[0-9]$")))
  expect_identical(rows(rulerx("m", "^[0-9]$", napass = FALSE)), 2L)
  expect_null(rows(rulemad("m", 10)))
  expect_identical(rows(rulemad("m", 10, napass = FALSE)), 2L)
  # nothing but missing values
  e <- data.frame(x = c(NA_real_, NA))
  expect_length(valrules(e, rulemad("x")), 0L)
  expect_identical(
    valrules(e, rulemad("x", napass = FALSE))[[1L]][["nfail"]],
    2L
  )
  expect_length(valrules(e, ruleincr("x")), 0L)
})

test_that("rulemad() is the median and the scaled median absolute deviation", {
  x <- c(10, 12, 11, 13, 12, 40, 9)
  med <- stats::median(x)
  dev <- 1.4826 * stats::median(abs(x - med))
  for (k in c(1, 2, 3, 20)) {
    got <- valrules(data.frame(x = x), rulemad("x", k))
    want <- which(abs(x - med) > k * dev)
    expect_identical(got[["x_within_mads"]][["rows"]], if (length(want)) want)
  }
  # no spread at the median: any other value fails
  z <- data.frame(x = c(5, 5, 5, 5, 6))
  expect_identical(valrules(z, rulemad("x"))[[1L]][["rows"]], 5L)
})

test_that("a rule about numbers on text says it could not be evaluated", {
  v <- valrules(d, list(rulebtw("div", 0, 1), ruleincr("div"), rulemad("div")))
  for (k in names(v)) {
    expect_match(
      v[[k]][["message"]],
      "could not be evaluated: the column is not numeric",
      fixed = TRUE
    )
  }
  expect_length(v, 3L)
  # dates and TRUE/FALSE are numbers here
  t <- data.frame(
    day = as.Date(c("2026-01-02", "2026-01-01")),
    ok = c(FALSE, TRUE)
  )
  expect_named(
    valrules(t, list(ruleincr("day"), ruleincr("ok"))),
    "day_increasing"
  )
  long <- data.frame(s = strrep("a", 5000), stringsAsFactors = FALSE)
  expect_match(
    valrules(long, rulerx("s", "a"))[[1L]][["message"]],
    "a value is longer than 4096 bytes",
    fixed = TRUE
  )
  expect_length(valrules(d[, 0], list(rulecmpl(), ruledstn())), 0L)
})

test_that("a rule written as data is the rule its function makes", {
  spec <- list(
    list(rule = "btw", col = "year", lo = 2000, hi = 2020),
    list(rule = "inset", col = "div", set = c("11", "14")),
    list(rule = "null", col = "m", sev = "fatal"),
    list(rule = "uniq", col = "year"),
    list(rule = "rx", col = "div", pat = "^[0-9]{2}$"),
    list(rule = "incr", col = "year", strict = TRUE),
    list(rule = "mad", col = "n"),
    list(rule = "mad", col = "n", n = 5),
    list(rule = "cmpl"),
    list(rule = "dstn", cols = "year"),
    list(rule = "colcnt", n = 3)
  )
  made <- list(
    rulebtw("year", 2000, 2020),
    ruleinset("div", c("11", "14")),
    rulenull("m", sev = "fatal"),
    ruleuniq("year"),
    rulerx("div", "^[0-9]{2}$"),
    ruleincr("year", strict = TRUE),
    rulemad("n"),
    rulemad("n", 5),
    rulecmpl(),
    ruledstn("year"),
    rulecolcnt(3)
  )
  expect_identical(valrules(d, spec), valrules(d, made))
  expect_identical(valrules(d, spec[[1L]]), valrules(d, made[[1L]]))
  # the stored form: kind, name, column, severity, message, what it was given
  expect_identical(
    unclass(rulebtw("year", 2000, 2020)),
    list(
      rule = "btw",
      name = "year_between",
      col = "year",
      sev = "warning",
      msg = "Column 'year' has values outside [2000, 2020]",
      lo = 2000,
      hi = 2020,
      napass = TRUE
    )
  )
  expect_identical(
    unclass(ruledstn()),
    list(
      rule = "dstn",
      name = "distinct_rows",
      sev = "warning",
      msg = "Some rows are duplicated"
    )
  )
  # what a written rule lacks is refused in the words of its function
  expect_error(
    valrules(d, list(list(rule = "btw", col = "year"))),
    "rulebtw: `lo` must be one number",
    fixed = TRUE
  )
  expect_error(
    valrules(d, list(list(rule = "null"))),
    "rulenull: `col` must be one non-empty string",
    fixed = TRUE
  )
  expect_error(
    valrules(d, list(list(rule = "inset", col = "div"))),
    "ruleinset: `set` must be a vector",
    fixed = TRUE
  )
})

test_that("the rules print", {
  expect_output(
    print(rulenull("m", sev = "fatal")),
    "rule m_not_null (fatal, on `m`): Column 'm' has missing values",
    fixed = TRUE
  )
  expect_output(print(rulecmpl()), "(warning, on the table)", fixed = TRUE)
  expect_invisible(print(rulecmpl()))
  v <- valrules(
    d,
    list(
      rulenull("m", sev = "fatal"),
      rulecolcnt(3),
      rulebtw("year", 2000, 2026)
    )
  )
  expect_output(print(v), "valrules: 3 rules, 2 failed (1 fatal)", fixed = TRUE)
  expect_output(
    print(v),
    "fatal  m_not_null: Column 'm' has missing values (1 of 4 rows); rows 2",
    fixed = TRUE
  )
  expect_output(
    print(v),
    "warning  col_count: Expected 3 columns",
    fixed = TRUE
  )
  expect_output(
    print(valrules(d, list())),
    "valrules: 0 rules, 0 failed",
    fixed = TRUE
  )
  big <- valrules(data.frame(x = 1:30), rulebtw("x", 100, 200))
  expect_output(
    print(big),
    "rows 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, ...",
    fixed = TRUE
  )
  expect_invisible(print(v))
})

test_that("the rule functions refuse bad input", {
  errs(
    "rul",
    list(
      list("`name` must be one non-empty string", "", function(v) TRUE),
      list("`name` must be one non-empty string", 1, function(v) TRUE),
      list("`expr` must be a function", "a", 1),
      list(
        "`col` must be one non-empty string",
        "a",
        function(v) TRUE,
        col = NA_character_
      ),
      list(
        "`msg` must be one non-empty string",
        "a",
        function(v) TRUE,
        msg = ""
      ),
      list("`sev` must be one of", "a", function(v) TRUE, sev = "bad")
    )
  )
  errs(
    "valrules",
    list(
      list("`data` must be a data frame", 1, list()),
      list("`rules` must be a rule or a list of rules", d, 3),
      list("`rules` must be a rule or a list of rules", d, list(1)),
      list(
        "`rules` must be a rule or a list of rules",
        d,
        list(list(rule = "zzz"))
      ),
      list(
        "`rules` must be a rule or a list of rules",
        d,
        list(list(rule = c("btw", "rx")))
      )
    )
  )
  errs(
    "ruleinset",
    list(
      list("`col` must be one non-empty string", "", "a"),
      list(
        "`set` must be a vector of at least one allowed value",
        "a",
        character(0)
      ),
      list(
        "`set` must be a vector of at least one allowed value",
        "a",
        list(1)
      ),
      list("`napass` must be TRUE or FALSE", "a", "b", napass = NA),
      list("`sev` must be one of", "a", "b", sev = "bad")
    )
  )
  errs(
    "rulebtw",
    list(
      list("`lo` must be one number", "a", "x", 1),
      list("`hi` must be one number", "a", 0, Inf),
      list("`lo` must not be above `hi`", "a", 2, 1),
      list("`napass` must be TRUE or FALSE", "a", 0, 1, napass = 2)
    )
  )
  errs(
    "rulenull",
    list(list("`col` must be one non-empty string", c("a", "b")))
  )
  errs("ruleuniq", list(list("`sev` must be one of", "a", sev = 1)))
  errs(
    "rulerx",
    list(
      list("`pat` must be one non-empty string", "a", ""),
      list("`pat` must be at most 256 bytes", "a", strrep("a", 257)),
      list("`pat` is not a valid regular expression", "a", "("),
      list("`napass` must be TRUE or FALSE", "a", "b", napass = NULL)
    )
  )
  errs(
    "ruleincr",
    list(list("`strict` must be TRUE or FALSE", "a", strict = 1))
  )
  errs(
    "rulemad",
    list(
      list("`n` must be one number", "a", n = "x"),
      list("`n` must be positive", "a", n = 0),
      list("`napass` must be TRUE or FALSE", "a", napass = NA)
    )
  )
  errs("rulecmpl", list(list("`sev` must be one of", sev = "bad")))
  errs(
    "ruledstn",
    list(
      list(
        "`cols` must be NULL or a character vector of column names",
        cols = 1
      ),
      list(
        "`cols` must be NULL or a character vector of column names",
        cols = character(0)
      ),
      list(
        "`cols` must be NULL or a character vector of column names",
        cols = c("a", "")
      )
    )
  )
  errs(
    "rulecolcnt",
    list(
      list("`n` must be one number", "x"),
      list("`n` must be a whole number of at least 0", 1.5),
      list("`n` must be a whole number of at least 0", -1)
    )
  )
})
