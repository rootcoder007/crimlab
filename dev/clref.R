# The reference values of the command line, the JSON entry clab_jrun()
# and the Python package for the functions they reach through clab_jrun():
# each case a call of the installed R package, written as a line of
# cli/test.tsv (function, arguments as JSON, the result as JSON, or ! and
# the message of a refusal). From the package root:
#   Rscript dev/clref.R [LIB]
# rewrites the lines of cli/test.tsv after the line that starts with
# "# made by dev/clref.R"; the lines above it are kept. LIB is a library
# holding crimlab.
#
# JSON as clab_jrun() reads and writes it, without a dependency: a data
# frame is an object of columns, a named vector an object of values, a
# list an object (an array when unnamed), an unnamed vector of one value
# outside a data frame that value, a matrix an array of its rows (an
# array of three dimensions nested, first index outermost); a factor is
# its labels, a date-time its seconds, a dist object its matrix. Doubles
# have 17 significant digits, NaN and infinities are strings, NA is null.
a <- commandArgs(TRUE)
if (length(a)) .libPaths(c(a[1L], .libPaths()))
suppressPackageStartupMessages(library(crimlab))

esc <- function(s) {
  s <- gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", enc2utf8(s)))
  # control characters as the command line writes them (\u000a), so that a
  # case stays on one line of test.tsv
  for (k in c(1:31, 127)) {
    s <- gsub(intToUtf8(k), sprintf("\\u%04x", k), s, fixed = TRUE)
  }
  paste0('"', s, '"')
}
one <- function(v) {
  if (is.character(v)) return(if (is.na(v)) "null" else esc(v))
  if (is.logical(v)) return(if (is.na(v)) "null" else if (v) "true" else "false")
  if (is.integer(v)) return(if (is.na(v)) "null" else sprintf("%d", v))
  if (is.nan(v)) return('"NaN"')
  if (is.na(v)) return("null")
  if (is.infinite(v)) return(if (v > 0) '"Inf"' else '"-Inf"')
  sprintf("%.17g", v)
}
arr <- function(x) paste0("[", paste(vapply(x, one, ""), collapse = ","), "]")
obj <- function(l) {
  if (!length(l)) return("{}")
  paste0("{", paste(sprintf("%s:%s", esc(names(l)), vapply(l, js, "")), collapse = ","), "}")
}
plain <- function(x) {
  # a clprep() result: its arguments, then its attributes
  if (inherits(x, "clpre")) {
    a <- attributes(x)
    b <- unclass(x)
    attributes(b) <- list(names = names(b))
    return(list(args = b, preset = a$preset, generic = a$generic, source = a$source,
                units = a$units, note = a$note, drop = a$drop))
  }
  if (is.factor(x)) {
    n <- names(x)
    x <- as.character(x)
    names(x) <- n
  }
  if (inherits(x, "POSIXct")) x <- as.numeric(x)
  if (inherits(x, "dist")) x <- as.matrix(x)
  x
}
js <- function(x) {
  x <- plain(x)
  if (is.null(x)) return("null")
  if (is.data.frame(x)) {
    cols <- vapply(x, function(v) arr(unname(as.vector(plain(v)))), "")
    return(paste0("{", paste(sprintf("%s:%s", esc(names(x)), cols), collapse = ","), "}"))
  }
  d <- dim(x)
  if (length(d) == 2L) {
    rows <- vapply(seq_len(d[1L]), function(i) arr(unname(x[i, ])), "")
    return(paste0("[", paste(rows, collapse = ","), "]"))
  }
  if (length(d) == 3L) {
    rows <- vapply(seq_len(d[1L]), function(i) js(matrix(x[i, , ], d[2L], d[3L])), "")
    return(paste0("[", paste(rows, collapse = ","), "]"))
  }
  if (is.list(x)) {
    if (!is.null(names(x))) return(obj(x))
    return(paste0("[", paste(vapply(x, js, ""), collapse = ","), "]"))
  }
  n <- names(x)
  x <- as.vector(x)
  names(x) <- n
  if (!is.null(names(x))) {
    return(paste0("{", paste(sprintf("%s:%s", esc(names(x)), vapply(unname(x), one, "")), collapse = ","), "}"))
  }
  if (length(x) == 1L) one(x) else arr(x)
}

cases <- character(0)
add <- function(fn, args, val) {
  cases[length(cases) + 1L] <<- paste(fn, obj(args), val, sep = "\t")
}
# A warning R gives is part of the result: an array "warnings" of its
# messages after the other elements of a list.
ok <- function(.f, ...) {
  args <- list(...)
  w <- character(0)
  r <- withCallingHandlers(do.call(.f, args), warning = function(m) {
    w <<- c(w, conditionMessage(m))
    invokeRestart("muffleWarning")
  })
  v <- js(r)
  if (length(w)) v <- sub("\\}$", paste0(',"warnings":', arr(w), "}"), v)
  add(.f, args, v)
}
bad <- function(.f, ...) {
  args <- list(...)
  m <- tryCatch({
    do.call(.f, args)
    NULL
  }, error = conditionMessage)
  if (is.null(m)) stop("no error from ", .f, call. = FALSE)
  add(.f, args, paste0("!", m))
}

# Fairness and the audits of risk tools
g8 <- rep(c("a", "b"), each = 4)
ok("fair", score = c(0.2, 0.7, 0.9, 0.4, 0.6, 0.1, 0.8, 0.3), y = c(0, 1, 1, 0, 1, 0, 0, 1), grp = g8)
ok("fair", score = c(0.2, 0.7, 0.9, 0.4, 0.6, 0.1, 0.8, 0.3), y = c(0, 1, 1, 0, 1, 0, 0, 1), grp = g8,
   cut = 0.35, base = "b", conf = 0.9)
bad("fair", score = c(0.2, 0.7), y = c(0, 2), grp = c("a", "b"))
ok("fairdf", pred = c(1, 0, 1, 1, 0, 0, 1, 0), grp = g8, y = c(1, 0, 1, 0, 0, 1, 1, 0))
bad("fairdf", pred = c(1, 0, 1, 1, 0, 0, 1, 0), grp = g8, y = c(1, 0, 1, 0, 0, 1, 1, 0), nboot = 1e9, seed = 1)
ok("fairdf", pred = c(1, 0, 1, 1, 0, 0, 1, 0), grp = g8, base = "b", nboot = 20, seed = 1)
bad("fairdf", pred = c(1, 0), grp = c("a", "a"))
g12 <- rep(c("a", "b"), 6)
ok("tmpaud", pred = c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1), grp = g12, period = rep(1:2, each = 6))
ok("tmpaud", pred = c(1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1), grp = g12, period = rep(1:2, 6),
   city = rep(c("x", "y"), each = 6), base = "b")
ok("lipaud", h = c(0.10, 0.20, 0.60, 0.65), d = dist(c(1, 2, 2.5, 6)), lip = 0.2)
ok("lipaud", h = c(0.1, 0.5, 0.4), pairs = rbind(c(1, 2), c(2, 3)), d = c(1, 2))
bad("lipaud", h = c(0.1, 0.5, 0.4))
ok("nsbnd", pobs = c(x = 0.1, y = 0.3, z = 0.5), a = 0.05, b = 0.1)
ok("nsbnd", pobs = c(0.2, 0.25), a = c(0, 0.1), b = 0.2)
bad("nsbnd", pobs = 0.3, a = 0.6, b = 0.5)
ok("rnkres", est = c(1, 2, 3.5), hw = 0.8)
ok("rnkres", est = c(a = 1, b = 2, c = 2.2, d = 5), hw = c(0.1, 0.5, 0.5, 1))
bad("rnkres", est = 1, hw = 1)
ok("nordet", x = rbind(c(0, 0), c(1, 1), c(3, 3)), off = data.frame(x = c(0, 2.5), y = c(0.5, 3)), radius = 1)
ok("nordet", x = rbind(c(0, 0), c(1, 1), c(3, 3)), off = rbind(c(0, 0.5)), radius = 2, p = 0.5, draw = TRUE, seed = 1)
ok("fairsim", n = 20, nareas = 4, seed = 1)
ok("fairsim", n = 15, groups = c("x", "y", "z"), props = c(1, 2, 1), nareas = 5, base = 0.4, bias = -0.2, seed = 2)
bad("fairsim", n = 10, groups = "A", seed = 1)
set.seed(3)
og <- rep(c("a", "b"), each = 20)
osev <- round(rnorm(40, ifelse(og == "a", 5, 4)), 2)
opri <- rpois(40, ifelse(og == "a", 2, 1))
olen <- round(2 + 0.8 * osev + 0.5 * opri + ifelse(og == "a", 1, 0) + rnorm(40), 2)
ok("oaxbd", y = olen, x = data.frame(sev = osev, pri = opri), grp = og, ref = "b")
ok("oaxbd", y = olen, x = osev, grp = og, conf = 0.9)
bad("oaxbd", y = olen, x = osev, grp = rep(c("a", "b", "c"), length.out = 40))

# Trace evidence
set.seed(1)
th <- matrix(round(rnorm(20, 10, 2), 3), 10, 2)
tx <- round(th[rep(1:10, each = 3), ] + matrix(rnorm(60, 0, 0.2), 30, 2), 3)
ok("tlvc", x = tx, item = rep(1:10, each = 3))
bg <- tlvc(tx, rep(1:10, each = 3))
ok("tlvlr", xa = th[1, ], xb = th[1, ] + 0.05, na = 3, nb = 3, bg = bg)
ok("tlvlr", xa = th[1, ], xb = th[2, ], na = 3, nb = 3, bg = bg)
bad("tlvlr", xa = 1, xb = 1, na = 3, nb = 3, bg = list(1))

# Aoristic weights, Benford, Scott-Knott, reporting panels
ok("aorwt", from = c(79200, 120600), to = c(93600, 124200), bw = 3600, nb = 24, origin = 0)
bad("aorwt", from = c(79200, 120600), to = c(93600, 124200), bw = 3600, nb = 1e9, origin = 0)
ok("aorwt", from = c(0, 3600, 7200, NA), to = c(1800, 10800, 7200, 9000), bw = 3600, nb = 4, origin = 0,
   each = TRUE)
bad("aorwt", from = c(0, 1), to = 1)
ok("bdig", x = c(123, 45.6, 0.29, -7, 0, NA, 9.99), nd = 1)
ok("bdig", x = c(123, 45.6, 0.29, -7, 1234), nd = 2, dp = 2, sgn = "both")
bad("bdig", x = 1, nd = 10)
ok("bnft", dig = c(1, 2, 1, 3, 1, 4, 1, 5, 9, 2, 1, 1, 6))
ok("bnft", dig = c(1, 2, 1, 3, 1, 4, 1, 5, 9, 2, 1, 1, 6), nsim = 49, seed = 1)
bad("bnft", dig = c(12, 1))
set.seed(1)
gx <- round(c(rnorm(5, 1.5180, 4e-5), rnorm(5, 1.5195, 4e-5)), 6)
ok("skp", x = gx, s2 = 4e-5^2, nsim = 99, seed = 1)
bad("skp", x = gx, s2 = 0, seed = 1)
ok("rptc", unit = rep(c("a", "b"), each = 3), year = rep(2001:2003, 2), rep = c(12, 12, 12, 5, 6, 7),
   y0 = 2001, y1 = 2003, minrep = 1)
bad("rptc", unit = "a", year = 2001.5, rep = 1, y0 = 2001, y1 = 2002)
dd <- c(1.2, 0.5, -0.3, 2, 1.1, 0.7, 0.9, -0.1)
ok("rswp", d = dd, gam = c(1, 1.5, 2))
ok("rswp", d = dd, alpha = 0.1)
bad("rswp", d = dd, gam = 0.5)
ok("rshl", d = dd, gam = c(1, 2))
ok("rshl", d = dd)
ok("lwss", x = 1:10, y = c(2, 4, 3, 5, 7, 6, 8, 9, 8, 10))
ok("lwss", x = c(3, 1, 2, 5, 4), y = c(1, 2, 2, 3, 5), f = 0.8, iter = 0, delta = 0)
bad("lwss", x = 1:3, y = 1:2)
set.seed(1)
cm <- cumsum(rnorm(330))
ok("chmb", a = round(cm[1:300] + rnorm(300, sd = 0.2), 4), b = round(cm[31:330] + rnorm(300, sd = 0.2), 4),
   wopt = 60, wval = 20)
bad("chmb", a = 1:10, b = 1:10, wopt = 5, wval = 6)
set.seed(1)
gm <- matrix(rpois(36, 2), 6)
gm[3:4, 3:4] <- 9
ok("gstr", cnt = gm)
ok("gstr", cnt = gm, ex = matrix(2, 6, 6), rad = 1, nsim = 19, seed = 1)
bad("gstr", cnt = gm, ex = matrix(1, 2, 2))
ok("sptst", a = c(5, 3, 8, 0, 2, 6), b = c(4, 6, 7, 1, 0, 9), nboot = 50, seed = 1)
ok("sptst", a = c(x = 5, y = 3, z = 8), b = c(z = 7, x = 4, y = 6), nboot = 20, seed = 2, alpha = 0.1)
bad("sptst", a = c(1, 2), b = 1, seed = 1)
tpar <- c(80, 50, 0.5, 0.05, 0.8, 0.9, 0.1, 0.2, 0.4, 0.6, 0.05, 0.1, 0.5, 0.7, 1, 1, 0.5)
ok("tfrp", par = tpar, ymax = 10, nsim = 100, seed = 1)
bad("tfrp", par = 1:3, seed = 1)
bad("tfrp", par = tpar, ymax = 1e9, nsim = 1, seed = 1)

# Trajectory models
set.seed(1)
zx <- cbind(1, (1:6 - 3.5) / 2.5)
zb <- cbind(c(0.5, -1), c(1.5, 0.5))
zy <- t(sapply(rep(1:2, each = 15), function(k) rpois(6, exp(zx %*% zb[, k]))))
ok("zgll", y = zy, x = zx, par = zb, prob = c(0.5, 0.5))
ok("zgem", y = zy, x = zx, par = zb + 0.2, prob = c(0.5, 0.5))
bad("zgll", y = zy, x = zx, par = zb, prob = c(0.5, 0.6))
set.seed(1)
mx <- cbind(1, 1:5)
my <- array(round(rnorm(12 * 5 * 2), 3), c(12, 5, 2))
my[1:6, , 1] <- my[1:6, , 1] + 3
mb <- array(c(0, 0, 0, 0, 3, 0, 0, 0), c(2, 2, 2))
ms <- array(diag(2), c(2, 2, 2))
ok("mvtll", y = my, x = mx, beta = mb, sig = ms, prob = c(0.5, 0.5))
ok("mvtem", y = my, x = mx, beta = mb, sig = ms, prob = c(0.5, 0.5))
bad("mvtll", y = my, x = mx, beta = mb, sig = ms, prob = c(0.2, 0.3, 0.5))

# Networks
nf <- c("a", "a", "b", "c", "d", "e", "e")
nt <- c("b", "c", "c", "d", "e", "f", "d")
ok("ntdeg", from = nf, to = nt)
ok("ntdeg", from = nf, to = nt, w = c(1, 2, 1, 3, 1, 2, 1), dir = TRUE, nodes = c(letters[1:6], "g"))
bad("ntdeg", from = nf, to = nt[1:3])
ok("ntbtw", from = nf, to = nt)
ok("ntbtw", from = nf, to = nt, dir = TRUE, norm = TRUE)
ok("ntcls", from = nf, to = nt)
ok("ntcls", from = nf, to = nt, dir = TRUE)
ok("ntevc", from = nf, to = nt)
ok("ntevc", from = nf, to = nt, w = c(1, 2, 1, 3, 1, 2, 1))
ok("ntcore", from = nf, to = nt)
ok("ntcmp", from = c(nf, "x"), to = c(nt, "y"), dir = TRUE)
ok("ntlvn", from = nf, to = nt, seed = 1)
ok("ntlvn", from = nf, to = nt, w = c(1, 2, 1, 3, 1, 2, 1), res = 0.5)
ok("ntmod", from = nf, to = nt, memb = c(a = 1, b = 1, c = 1, d = 2, e = 2, f = 2))
ok("ntmod", from = nf, to = nt, memb = c(1, 1, 1, 2, 2, 2), res = 0.5)
bad("ntmod", from = nf, to = nt, memb = c(1, 2))
ok("ntbrk", from = nf, to = nt, grp = c(a = "x", b = "x", c = "y", d = "y", e = "z", f = "z"))
ok("ntbrk", from = nf, to = nt, grp = c("x", "x", "y", "y", "z", "z"), dir = TRUE)
ok("ntchg", from = nf, to = nt)
bad("ntchg", from = c(nf, "x"), to = c(nt, "y"))
ok("coofg", inc = c(1, 1, 2, 2, 2, 3, 3), off = c("p", "q", "p", "r", "s", "q", "r"))
bad("coofg", inc = 1:2, off = "p")

# Checks
ok("relchk", x = data.frame(yes = c(30, 20), no = c(70, 80)), value = 1.71)
ok("relchk", x = list(yes = c(m = 30, f = 20), no = c(m = 70, f = 80)), value = 0.1, type = "rd", group = "m")
ok("relchk", x = rbind(c(30, 70), c(20, 80)), value = 1.5, type = "rr", tol = 0.01)
bad("relchk", x = rbind(c(30, 70), c(20, 80)), value = 1.5, type = "hr")
ok("catchk", x = data.frame(a = c("x", "y", "x", "Y", NA), b = c(0, 1, 1, 0, 1), c = c(1, 2, 3, 4, 5)))
ok("catchk", x = c("u", "v", "v", "w", "w", "w"), ref = "v", nmin = 2)
ok("catchk", x = list(s = c(TRUE, FALSE, TRUE), t = c(1.5, 2, 3)), ref = c(s = "TRUE"))
bad("catchk", x = list(a = 1, a = 2))

# Geography and centrography
ok("hvdist", from = c(-79.38, 43.65), to = rbind(c(-73.57, 45.50), c(-123.12, 49.28)))
bad("hvdist", from = c(-79.38, 43.65), to = c(200, 0))
ll <- data.frame(lon = c(-79.40, -79.38, -79.36, -79.39), lat = c(43.64, 43.66, 43.65, 43.67))
ok("eqproj", x = ll)
ok("eqproj", x = ll$lon, y = ll$lat, lon0 = -79.4, lat0 = 43.6)
ok("eqinv", x = c(0, 1000, -500), y = c(0, 2000, 300), lon0 = -79.4, lat0 = 43.6)
pxy <- data.frame(x = c(0, 1, 2, 3, 4, 2.5), y = c(0, 2, 1, 3, 0.5, 4))
ok("cgstat", x = pxy)
ok("cgstat", x = pxy$x, y = pxy$y, w = c(1, 2, 1, 1, 3, 1), k = 0)
ok("cgstat", x = ll, lonlat = TRUE)
bad("cgstat", x = ll)
ok("cgmedn", x = pxy)
ok("webpt", x = pxy)
ok("webpt", x = pxy, w = c(1, 2, 1, 1, 3, 1), start = c(1, 1), tol = 1e-8, maxit = 50)
ok("sdell", x = pxy)
ok("sdell", x = pxy, k = 2, c = 2)
ok("cvhull", x = pxy)

# Point patterns, hot spots
set.seed(2)
pp <- data.frame(x = round(runif(30, 0, 10), 3), y = round(runif(30, 0, 10), 3))
w10 <- c(0, 10, 0, 10)
ok("nnidx", x = pp, win = w10)
ok("nnidx", x = pp, win = w10, kmax = 2, nsim = 19, seed = 1, dist = TRUE, metric = "manhattan")
bad("nnidx", x = pp, win = c(0, 1, 0, 1))
ok("ripk", x = pp, win = w10, r = c(0.5, 1, 2))
ok("ripk", x = pp, win = w10, r = c(1, 2), corr = "border", metric = "manhattan")
bad("ripk", x = pp, win = w10, r = c(2, 1))
ok("ripenv", x = pp, win = w10, r = c(0.5, 1, 2), nsim = 19, seed = 1)
ok("ripk", x = pp, win = w10, r = c(5e-324, 1e9))
ok("ripenv", x = pp, win = w10, r = c(5e-324, 1e9), nsim = 19, seed = 1)
ok("fzmode", x = pp, r = 2)
ok("fzmode", x = pp, r = 2, count = rep(1:3, 10), metric = "manhattan")
ok("kmpart", x = pp, k = 3, seed = 1)
ok("kmpart", x = pp, k = 2, init = rbind(c(2, 2), c(8, 8)), upd = "weber")
ok("kmpart", x = pp, k = 2, init = "cells", cell = c(5, 5), count = rep(c(1, 2), 15))
set.seed(4)
hp <- data.frame(x = round(c(rnorm(25, 2, 0.3), runif(25, 0, 10)), 3), y = round(c(rnorm(25, 2, 0.3), runif(25, 0, 10)), 3))
ok("nnhclu", x = hp, area = 100, minpts = 5)
ok("nnhclu", x = hp, fixd = 0.5, minpts = 4, maxord = 2, count = rep(1:2, 25))
bad("nnhclu", x = hp)
ok("stacel", x = hp, win = w10, r = 1, top = 5, minpts = 4)
ok("stacel", x = hp, win = w10, nx = 10, ny = 8, r = 1.5, top = 3, minpts = 3, count = rep(1:2, 25))
ok("rnnhcl", x = hp, base = pp, h = 2, minpts = 5)
ok("rnnhcl", x = hp, base = pp, bw = rep(2, 30), h = 1, p = 0.3, minpts = 4, k = 0)

# Space and time
set.seed(5)
st <- data.frame(x = round(runif(20, 0, 5), 3), y = round(runif(20, 0, 5), 3))
stt <- round(sort(runif(20, 0, 30)), 2)
ok("mantst", x = st, t = stt, nsim = 49, seed = 1)
ok("mantst", x = st, t = stt, cs = 1, ct = 1, nsim = 0, metric = "manhattan")
bad("mantst", x = st[1:2, ], t = stt[1:2])
ok("stmavg", x = st, t = stt)
ok("stmavg", x = st, t = stt, span = 3)
bad("stmavg", x = st, t = stt, span = 4)

# Journey to crime
ok("jtdecay", d = c(0, 0.5, 1, 2, 5), type = "negexp", par = c(1, 0.5))
ok("jtdecay", d = c(0, 0.5, 1, 2, 5), type = "table", par = 1, tab = c(1, 0.8, 0.5, 0.2))
bad("jtdecay", d = 1, type = "linear", par = 1)
ok("jtsurf", x = c(2, 4, 3), y = c(3, 1, 4), win = c(0, 6, 0, 6), nx = 4, ny = 3, type = "negexp", par = c(1, 0.4))
ok("jtsurf", x = c(2, 4, 3), y = c(3, 1, 4), win = c(0, 6, 0, 6), nx = 3, ny = 3, type = "linear",
   par = c(1, 0.1), metric = "manhattan", mode = "lik")
ok("jtcal", d = c(0.5, 1, 1.5, 2, 3, 4), t = c(0, 1, 2, 3), h = 0.8)
ok("jtcal", d = c(0.5, 1, 1.5, 2, 3, 4), t = c(-1, 1, 2), h = 0.5, w = c(1, 2, 1, 1, 1, 1), kern = "quartic",
   refl = FALSE)
js1 <- jtsurf(c(2, 4, 3), c(3, 1, 4), c(0, 6, 0, 6), 4, 3, "negexp", c(1, 0.4))
ok("jtacc", prob = js1$prob, win = c(0, 6, 0, 6), base = c(3, 2))
bad("jtacc", prob = js1$prob, win = c(0, 6, 0, 6), base = c(9, 2))
ok("bjtc", pj = c(0.2, 0.5, 0.3), od = rbind(c(2, 1), c(1, 3), c(0, 2)), dest = c(1, 2))
bad("bjtc", pj = c(0.2, 0.5, 0.3), od = rbind(c(2, 1)), dest = 1)

# Spatial weights and autocorrelation
ok("spwt", x = pp[1:8, ], k = 3)
ok("spwt", x = pp[1:8, ], type = "band", d = 4, idw = 1, style = "raw")
ok("spwt", type = "contig", nb = list(2, c(1, 3), c(2, 4, 5), 3, 3))
bad("spwt", type = "knn")
wk <- spwt(pp, k = 4)
xv <- round(pp$x + pp$y / 2, 3)
ok("gearc", x = xv, w = wk)
ok("gearc", x = xv, w = wk, nsim = 19, seed = 1)
ok("gorg", x = xv, w = wk, nsim = 19, seed = 1)
ok("lmorn", x = xv, w = wk)
ok("lmorn", x = xv, w = wk, nsim = 19, seed = 1)
ok("locgi", x = xv, w = wk)
ok("locgi", x = xv, w = wk, star = FALSE)
bad("locgi", x = xv[1:5], w = wk)
ok("dirmn", angle = c(0.1, 0.4, -0.2, 1))
ok("dirmn", x = pxy[-1, ], origin = c(0, 0), w = c(1, 2, 1, 1, 1))
bad("dirmn", x = pxy, angle = 1)
ok("hdbang", x = pp, v = xv, nn = 5)
ok("hdbang", x = pp, v = xv, w = rep(c(1, 2), 15), nn = 4, ang = 120, ntrip = 2, iter = 2)
set.seed(6)
cw <- data.frame(x = round(cumsum(rnorm(12)), 3), y = round(cumsum(rnorm(12)), 3))
ok("crwalk", x = cw, h = 3, nsim = 9, seed = 1)
bad("crwalk", x = cw[1:2, ], seed = 1)

# Forensic toxicology
ok("txcal", conc = c(1, 2, 5, 10, 20), resp = c(0.11, 0.2, 0.52, 0.98, 2.05), unk = c(0.05, 0.4, 1.5, NA))
ok("txcal", conc = c(1, 2, 5, 10, 20), resp = c(0.11, 0.2, 0.52, 0.98, 2.05), w = "none", k = c(3, 9))
bad("txcal", conc = c(1, 2), resp = c(1, 2))
ok("lodimp", x = c(0.4, 0.05, NA, 1.2, 0.08), lod = 0.1)
ok("lodimp", x = c(0.4, 0.05, NA, 1.2, 0.08), lod = c(0.5, 0.01, 0.1, 0.1, 0.1), meth = "sqrt2")
ok("pmrrat", cen = c(0.5, 1, 2, 0.3), per = c(0.4, 0.6, 0.5, 0.3))
ok("pmrrat", cen = c(0.5, 1, 2), per = c(0.4, 0.6, 0.5), cut = c(1.5, 3))
bad("pmrrat", cen = 1, per = 1, cut = c(2, 1))
ok("etohlr", eth = c(0.5, 0.1, 0.8), prop = c(0, 0.02, 0), but = 0, etg = c(1, NA, 0), thr = 0.01)
ok("etohlr", eth = c(0.5, 0.1))
ok("mtxrel", base = c(blood = 0.9, vitreous = 0.8, urine = 0.6), dil = c(TRUE, FALSE, FALSE), days = 3, stage = 2)
bad("mtxrel", base = c(0.9, 0.8), dil = TRUE)

# Deterrence and severity
ok("detrsp", x = c(0, 1, 2, 3), ben = c(0, 1.5, 2.5, 3), sanc = c(0, 1, 3, 6), p = c(0.1, 0.5, 0.9))
bad("detrsp", x = c(0, 1, 1), ben = c(0, 1, 2), sanc = c(0, 1, 2), p = 0.5)
dp <- c(0.1, 0.2, 0.3, 0.4, 0.2, 0.3)
ds <- c(1, 2, 1, 3, 2, 2)
dc <- c(0.5, 1, 1.5, 0.5, 2, 1)
ok("dtrid", p = dp, s = ds, c = dc)
ok("dtrid", p = dp, s = ds, c = dc, y = c(3, 2.5, 2.8, 1.5, 2.1, 2), log = TRUE, rho = 0.1)
ok("inspg", a = 0.5, x = c(0.5, 0.3, 0.2), steps = 5, dt = 0.1)
bad("inspg", a = 0.5, x = c(0.5, 0.3, 0.3))
bad("inspg", a = "a", x = c(0.5, 0.3, 0.2))
bad("dtrid", p = dp, s = ds, c = dc, rho = "a")
ok("insphs", temp = c(0.1, 0.5), cost = c(0.2, 0.8, 1), steps = 20)
ok("colarr", a = c(0.2, 0.5), e = 0.3, bg = c(0, 0.5))
bad("colarr", a = c(0.2, 0.5), e = c(0.3, 0.2, 0.1), bg = 0)
ok("crsev", x = c(2, 1, 3, 1), off = c("theft", "assault", "theft", "robbery"),
   w = c(assault = 5, robbery = 4, theft = 1))
ok("crsev", x = c(2, 1, 3, 1, 2, 2), off = c("t", "a", "t", "a", "t", "a"), w = c(1, 5),
   unit = rep(c("u", "v"), each = 3), time = c(1, 1, 2, 1, 2, 2), pop = c(100, 100, 200, 50, 60, 60),
   base = 1, sent = c(a = 2, t = 1))
bad("crsev", x = 1, off = "a", w = c(b = 1))
ok("zcomp", x = c(3, 1, 4, 1, 5, 9, 2, 6), unit = rep(c("p", "q", "r", "s"), 2), cat = rep(c("x", "y"), each = 4))
ok("zcomp", x = c(3, 1, 4, 1, 5, 9, 2, 6), unit = rep(c("p", "q", "r", "s"), 2), cat = rep(c("x", "y"), each = 4),
   w = c(y = 2, x = 1))

# Concentration and the dark figure
cx <- c(0, 0, 1, 0, 3, 7, 0, 1, 2, 0, 12, 1)
ok("cconc", x = cx, nsim = 0)
bad("cconc", x = cx, nsim = 1e9, seed = 1)
ok("cconc", x = cx, w = rep(c(1, 2), 6), by = rep(c("n", "s"), each = 6), fracs = c(0.1, 0.5), pop = 1:12 * 100, nsim = 0)
ok("cconc", x = c("a", "b", "a", "c", "a", "b", "d"), by = c(1, 1, 1, 2, 2, 2, 2), nsim = 19, seed = 1)
bad("cconc", x = c(0, 0))
set.seed(7)
hx <- round(1 / runif(80)^(1 / 1.5))
ok("hilltl", x = hx)
ok("hilltl", x = hx, xmin = 3, discrete = FALSE)
ok("hilltl", x = c(1, 2, 3), mintl = 5)
bad("hilltl", x = c(0, 0))
ok("crcap", x = c(n1 = 120, n2 = 90, m = 30))
ok("crcap", x = c("100" = 40, "010" = 30, "001" = 25, "110" = 10, "101" = 8, "011" = 6, "111" = 3),
   kappa = 1.5, int = c(1, 2))
bad("crcap", x = c(1, 2))
ok("dkfig", v = c(0.1, 0.2), r = 0.05, amax = 0.1, bmax = 0.05, thr = 0.3)
ok("dkfig", va = 0.2, ra = 0.05, rb = 0.02, lam = 0.5, k = c(1, 1, 2, 3))
bad("dkfig", amax = 0.1)

# Crosswalks and taphonomy
ok("bndxw", from = c("a", "a", "b", "c"), to = c("X", "Y", "Y", "Y"), x = c(10, 20, 5, 1), unit = c("a", "b", "c", "z"),
   w = c(0.6, 0.4, 1, 1))
ok("bndxw", from = c("a", "b"), to = c("X", "X"), x = c(10, 20), unit = c("a", "b"))
bad("bndxw", from = "a", to = "X", x = 1, unit = c("a", "b"))
set.seed(8)
sx <- data.frame(age = round(runif(40, 15, 60)), male = rbinom(40, 1, 0.5))
sy <- rbinom(40, 1, plogis(-2 + 0.03 * sx$age + 0.5 * sx$male))
ok("synexp", y = sy, x = sx, ax = data.frame(age = c(20, 40), male = c(0.5, 0.5)), pop = c(1000, 2000))
bad("synexp", y = rep(1, 40), x = sx, ax = sx[1:2, ], pop = c(1, 1))
cp <- rbind(c(0.2, 0.3, 0.5), c(0.1, 0.1, 0.8), c(0, 0.5, 0.5))
ok("clrilr", x = cp)
ok("clrilr", x = data.frame(a = c(1, 2), b = c(3, 4), c = c(5, 6)), type = "ilr", ps = 0)
bad("clrilr", x = rbind(c(0, 1)), ps = 0)
ok("preslr", x = c(10, 12, 15), m1 = 11, s1 = 2, m2 = 14, s2 = c(2, 2, 3))
ok("tpdmc", pres = 0.2, decay = 0.4, mum = 0.3)
ok("tpdmc", states = c("fresh", "bloat"), start = "bloat", steps = 10, seed = 1)
bad("tpdmc", states = c("a", "a"))
bad("tpdmc", steps = 2147483647, seed = 1)
set.seed(9)
tby <- round(rnorm(30, 10), 3)
tbx <- data.frame(temp = round(runif(30, 5, 30), 2))
ok("tpbhm", y = tby, x = tbx)
ok("tpbhm", y = tby, x = tbx, grp = rep(c("p", "q", "r"), 10), pm = c(0, 0.1), psd = 5, level = 0.9)
bad("tpbhm", y = tby, x = tbx[1:5, , drop = FALSE])
ok("dirsim", n = 3, alpha = c(30, 15, 5), seed = 1)
ok("dirsim", n = 2, alpha = c(85, 5, 2, 0.5), total = 1e6, seed = 7)
bad("dirsim", n = 2, alpha = c(1, 0), seed = 1)

# Learners
set.seed(10)
lx <- matrix(round(rnorm(60), 3), 20, 3)
ly <- round(lx[, 1] - 0.5 * lx[, 2] + rnorm(20, sd = 0.3), 3)
lb <- rbinom(20, 1, plogis(0.8 * lx[, 1]))
ok("enpath", x = lx, y = ly, nlam = 5)
bad("enpath", x = lx, y = ly, nlam = 1e9)
ok("enpath", x = lx, y = lb, fam = "binomial", alpha = 0.5, lambda = c(0.2, 0.1, 0.05), fold = rep(1:4, 5))
bad("enpath", x = lx, y = rep(1, 20))
ok("lgtreg", x = lx, y = lb)
bad("lgtreg", x = lx[1:3, ], y = c(0, 1, 0))
ok("rfor", x = lx, y = ly, ntree = 20, seed = 2)
ok("rfor", x = lx, y = ly, newx = lx[1:3, ], ntree = 10, mtry = 2, minnode = 3, frac = 0.8, repl = FALSE, hon = TRUE)
bad("rfor", x = lx, y = ly, mtry = 9)

# Profiling, risk terrain and the life course
ok("crchyp", x = c(0, 4, 2, 1), y = c(0, 0, 3, 1), home = c(2, 1))
ok("crchyp", x = c(0, 4, 2, 1), y = c(0, 0, 3, 1), home = rbind(c(2, 1), c(9, 9)))
bad("crchyp", x = 1, y = 1)
ok("hitscr", s = c(0.1, 0.5, 0.5, 0.9, 0.2), home = c(2, 4))
bad("hitscr", s = c(0.1, 0.5), home = 3)
set.seed(11)
ry <- rpois(25, 2)
ok("rtmfit", y = ry, lay = data.frame(bars = round(runif(25), 3), parks = rbinom(25, 1, 0.3)))
rc <- as.matrix(expand.grid(x = 1:5, y = 1:5))
ok("rtmfit", y = ry, lay = list(bars = rbind(c(2, 2), c(4, 4))), cell = rc, rad = 1.5, kind = "dens", level = 0.9)
bad("rtmfit", y = rep(0, 3), lay = 1:3)
ok("agemix", curves = data.frame(early = c(0.1, 0.3, 0.2, 0.1), late = c(0, 0.1, 0.2, 0.3)), shares = c(0.4, 0.6),
   ages = c(15, 20, 25, 30))
ok("agemix", panel = cbind(c(1, 0, 2, 3, 0, 1), c(0, 0, 1, 2, 1, 1)), nperm = 99, seed = 1)
bad("agemix", shares = 1)
ok("cumrsk", age = c(10, 15, 20, 25), ev = c(5, 20, 15, 5), pop = c(1000, 1000, 950, 900))
ok("cumrsk", age = c(10, 15, 20, 25), ev = c(5, 20, 15, 5), pop = c(1000, 1000, 950, 900), atrisk = c(990, 970, 930, 900),
   dth = c(0.001, 0.002, 0.002, 0.003), level = 0.9)
bad("cumrsk", age = c(2, 1), ev = c(1, 1), pop = c(10, 10))
ok("jpipe", num = c(5, 6, 7, 8, 6, 5), den = c(10, 12, 11, 15, 12, 10), win = 3)
ok("jpipe", num = c(5, 6, 7, 8, 6, 5), den = c(10, 12, 11, 15, 12, 10), win = 2, lag = 1, grp = c("a", "b", "a", "b", "a", "b"))

# Covariate balance and the E-value
set.seed(12)
bx <- data.frame(age = round(rnorm(30, 40, 10)), prior = rpois(30, 2))
bd <- rbinom(30, 1, plogis(-2 + 0.05 * bx$age))
ok("smdbal", x = bx, d = bd)
ok("smdbal", x = unname(as.matrix(bx)), d = bd, w = rep(c(1, 2), 15), ps = rep(c(0.3, 0.6), 15), thr = 0.2)
bad("smdbal", x = bx, d = rep(0, 30))
ok("evlx", est = 1.8, lo = 1.2, hi = 2.6)
ok("evlx", est = 0.4, type = "smd", sd = 0.1)
ok("evlx", est = 2.2, type = "or", rare = TRUE)
ok("evlx", est = 3.9, lo = 1.8, hi = 8.7, true = 2)
bad("evlx", est = 2, true = 0)
bad("evlx", est = 1.8, lo = 2, hi = 3)

# Courts
set.seed(13)
jz <- rbinom(40, 1, 0.5)
jd <- rbinom(40, 1, 0.3 + 0.3 * jz)
jy <- round(rnorm(40, 1 + jd), 3)
ok("jdgiv", z = jz, d = jd, y = jy)
ok("jdgiv", z = jz, d = jd, y = jy, w = rep(c(1, 2), 20), dsh = c(0, 0.1), deff = c(0, 1))
ok("jdgiv", d = cbind(c(0, 1, 1, 0, 1), c(1, 1, 1, 0, 0)), y = data.frame(a = c(1, 2, 3, 4, 5), b = c(2, 2, 4, 4, 5)))
bad("jdgiv", z = rep(1, 4), d = c(0, 1, 0, 1), y = 1:4)
ok("mmshr", sent = c(12, 24, 36, 12, 60), min = 24)
ok("mmshr", sent = c(12, 24, 36, 12, 60), min = c(12, 24, 24, 12, 48), grp = c("a", "b", "a", "b", "a"))
bad("mmshr", sent = c(1, 2, 3), min = c(1, 2))

# Self-exciting processes and change points
ht <- c(0.5, 1.2, 1.3, 2.8, 3.1, 3.15, 5.6, 7.2, 7.25, 7.3, 9.1)
ok("hwkres", times = ht, hor = 10, par = c(0, 0.4, 2))
ok("hwkres", times = ht, hor = 10, par = c(0, 0.2, 0.1, 0.1, 0.3, 1.5, 2), kernel = "weibull", base = "cycle", per = 5)
bad("hwkres", times = rev(ht), hor = 10, par = c(0, 0.4, 2))
ok("hwksim", par = c(0, 0.4, 2), hor = 20, seed = 1)
ok("hwksim", par = c(0, 0.4, 2), hor = 10, win = c(0, 1, 0, 1), sigma = 0.1, seed = 2)
bad("hwksim", par = c(0, 2, 2), hor = 100, nmax = 5, seed = 1)
bad("hwksim", par = c(0, 0.4, 2), hor = 20, nmax = 1e9, seed = 1)
sy2 <- c(3, 4, 3, 5, 4, 9, 10, 8, 11, 9)
ok("stpcp", y = sy2, nperm = 99, seed = 1)
ok("stpcp", y = sy2, x = 2001:2010, minseg = 3, nperm = 49, seed = 2)
bad("stpcp", y = sy2, x = 1:3, seed = 1)

# Counts, rates and tables
# counts spread well beyond Poisson, so the negative binomial's dispersion is identified and its fit
# stops at the same point on every platform
cr <- data.frame(n = c(12, 60, 3, 25, 40, 141, 2, 22, 11, 95, 7, 8), female = rep(c("no", "yes"), 6),
                 year = rep(c("2021", "2022", "2023"), each = 4), inst = rep(c("A", "B", "C", "D"), each = 2, length.out = 12))
ok("cntreg", y = cr$n, grp = cr$female, adj = cr["year"], clus = cr$inst)
ok("cntreg", y = c(5, 12), grp = c("a", "b"), expo = c(1000, 1500))
ok("cntreg", y = cr$n, grp = cr$female, adj = cbind(seq(0.5, 6, 0.5)), clus = cr$inst, gee = TRUE, cadj = FALSE, ref = "yes",
   level = 0.9)
bad("cntreg", y = c(0, 0, 0), grp = c("a", "b", "a"))
ok("grprr", x = c(white = 120, black = 95, other = 30), n = c(white = 60000, black = 21000, other = 15000), per = 1000)
ok("grprr", x = c(1, 0, 0, 1, 1, 0, 1, 1), grp = c(1, 1, 1, 1, 2, 2, 2, 2), type = "prop", ref = 2, level = 0.9)
bad("grprr", x = c(1, 2), type = "rate")
da <- data.frame(age = rep(c("18-29", "30-49", "50+"), 2), grp = rep(c("A", "B"), each = 3), x = c(30, 22, 5, 41, 30, 12),
                 n = c(5000, 9000, 6000, 3000, 8000, 7000))
ok("dastd", x = da$x, n = da$n, strat = da$age, grp = da$grp, std = c(30, 45, 25), per = 1e4)
ok("dastd", x = da$x, n = da$n, strat = da$age, grp = da$grp, std = c("50+" = 1, "18-29" = 2, "30-49" = 3), ref = "B")
bad("dastd", x = da$x, n = da$n, strat = da$age, grp = da$grp, std = c(1, 2))
si <- data.frame(area = rep(c("N", "S", "E"), each = 2), age = rep(c("young", "old"), 3), x = c(12, 3, 30, 9, 4, 1),
                 n = c(4000, 2000, 8000, 3000, 1500, 900))
ok("sirex", x = si$x, n = si$n, strat = si$age, grp = si$area)
ok("sirex", x = si$x, n = si$n, strat = si$age, grp = si$area, rates = c(old = 0.003, young = 0.002), level = 0.9)
bad("sirex", x = si$x, n = si$n, strat = si$age, grp = si$area, rates = c(old = 0, young = 0))
ok("ctabx", row = c("none", "none", "knife", "gun", "none", "knife", "gun", "none"),
   col = c("hands", "hands", "CEW", "firearm", "CEW", "hands", "firearm", "hands"))
ok("ctabx", row = c("N", "N", "E", "W", "E"), col = c("N", "E", "E", "W", "E"))
ok("ctabx", row = c("A", "A", "B", "B"), col = c("F", "M", "F", "M"), w = c(12, 40, 9, 51), exact = TRUE)
ok("ctabx", row = rbind(c(10, 4), c(3, 12)), correct = TRUE, level = 0.9)
ok("ctabx", row = c("a", "a", "b", "b", "b"), col = c("x", "y", "x", "x", "y"), w = c(3, 1, 4, 1, 5), agg = "mean")
bad("ctabx", row = c("a", "b"), col = c("x", "y"), agg = "max")
bad("ctabx", row = matrix(0, 2, 2))
bad("wrtxtfb", text = "a", path = strrep("x", 5000))

# Audits of risk scores
set.seed(31)
ca <- data.frame(area = rep(paste0("a", 1:6), each = 8), risk = round(runif(48), 3), group = rep(c("g", "h"), 24))
ca$outcome <- rbinom(48, 1, ca$risk)
ok("calaud", area = ca$area, risk = ca$risk, outcome = ca$outcome, group = ca$group, nperm = 49, seed = 1)
ok("calaud", area = ca, nperm = 0)
ok("calaud", area = ca$area, risk = ca$risk, outcome = ca$outcome, group = ca$group, nperm = 0,
   population = setNames(c(100, 200, 150, 120, 90, 300), paste0("a", 1:6)))
bad("calaud", area = c("a", "b"), risk = c(0.1, 0.2), outcome = c(0, 1), group = c("g", "h"), nperm = 0)
set.seed(32)
sg <- rep(c("x", "y", "z"), c(8, 8, 6))
ss <- round(rnorm(22, ifelse(sg == "x", 0, ifelse(sg == "y", 0.5, 1))), 3)
ok("scrdsp", score = ss, group = sg)
ok("scrdsp", score = ss, group = sg, ref = "z", conf = 0.9)
bad("scrdsp", score = c(1, 2), group = c("a", "a"))

# Courts: judges' leniency
set.seed(33)
jj <- sample(paste0("J", 1:4), 80, TRUE)
jd <- rbinom(80, 1, c(J1 = 0.2, J2 = 0.4, J3 = 0.5, J4 = 0.7)[jj])
jy <- round(rnorm(80, jd), 3)
ok("jmono", judge = jj, d = jd)
ok("jmono", judge = jj, d = jd, y = jy, grp = rep(c("adult", "youth"), 40), minn = 5, adj = "BH")
bad("jmono", judge = c("a", "b"), d = c(0, 1), w = c(0, 0))

# Long stays
set.seed(34)
md <- data.frame(id = sample(1:15, 30, TRUE), yr = sample(2021:2022, 30, TRUE), days = rgeom(30, 0.08) + 1)
ok("mndla", days = md$days, id = md$id, by = md$yr)
ok("mndla", days = c(5, 20, 30), hours = c(1, 1.5, 3), miss = c(100, 100, 40))
ok("mndla", days = c(3, 8, 20, 25), alert = data.frame(a = c(TRUE, FALSE, TRUE, TRUE), b = c("yes", "no", "no", "yes")), amin = 1,
   elig = c(1, 1, 0, 1), thr = 10, level = 0.9)
bad("mndla", days = c(3, 8), hours = c(1, 2))

# Published figures
yd <- data.frame(year = rep(2019:2023, 2), svc = rep(c("A", "B"), each = 5), n = c(120, 131, 98, 105, 140, 30, 28, 41, 25, 33))
ok("yoypx", x = yd$n, prd = yd$year, by = yd$svc, rnd = 5)
ok("yoypx", x = c(10.5, 12, 9, 15), prd = 1:4, units = "percent", minbase = 10)
ok("yoypx", x = c(4, 0, 6, 8), prd = c(1, 2, 3, 4), lag = 2, adj = "holm", alpha = 0.1)
bad("yoypx", x = c(3e50, 5e50), prd = 1:2)
bad("opteffort", poc = c(0.5, 0.5), area = c(1, 1), W = 5e-324, effort = 1)
bad("yoypx", x = c(1, 2), prd = c(1, 2), units = "continuous", rnd = 5)
ok("bandx", x = c("1-5", "6-15", "16-30", "31-60", "Over 60"), cnts = c(4, 2, 1, 1, 2))
ok("bandx", x = c("Under 18", "18 to 24", "25-34", "35 and over"), rule = "geometric", cap = 80, fac = 1.5, floor = 12)
ok("bandx", x = c("0-4", "5-9", "lots", "10+"), cnts = c(3, 2, 0, 1), clsup = FALSE)
bad("bandx", x = c("1-5", "lots"), cnts = c(1, 1))
ok("pbbnd", x = c("1,240", "35", "<5", "x"), rnd = 5, suplim = 10)
ok("pbbnd", x = c(120, 15), rnd = 5, prev = c(100, 20), pop = 5e4)
ok("pbbnd", x = c(120, 15), rnd = 5, prev = c("100", "x"), pop = c(5e4, 6e4), prepop = 4e4, suplim = 5, per = 100)
bad("pbbnd", x = c("12", "abc"))

# Stocks and flows
set.seed(35)
ls <- round(rlnorm(40, 2.3, 1.1)) + 1
ok("lostay", days = ls, breaks = c(0, 5, 15, 31, 61, Inf), thr = 15)
ok("lostay", days = c(3, 8, 20, 25, 40, 2), w = c(2, 1, 1, 3, 1, 1), flag = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE),
   breaks = c(0, 10, 50), thr = c(5, 30), level = 0.9)
bad("lostay", days = c(1, 2), w = c(0.5, 1))
ok("dktflw", recv = c(120, 130, 125, 140), disp = c(110, 135, 120, 150), pend = c(400, 395, 400, 390), win = 2, span = 91)
ok("dktflw", mix = rbind(c(50, 30, 20), c(55, 28, 17), c(60, 25, 15)))
ok("dktflw", recv = c(120, 130, 125), disp = c(110, 135, 120), pend = c(400, 395, 400), prd = c("q1", "q2", "q3"),
   span = c(90, 91, 92), mtime = c(300, 290, 310), shares = rbind(c(0.5, 0.5), c(0.6, 0.4), c(0.4, 0.6)),
   floors = c(0, 90))
bad("dktflw", recv = c(1, 2))
ok("stkflw", people = c(5200, 4800, 4300), days = c(310000, 318000, 330000), prd = 2021:2023)
ok("stkflw", adp = 900, alos = 61)
ok("stkflw", people = c(500, 520), alos = c(30, 28), from = c("2021-01-01", "2022-01-01"), to = c("2021-12-31", "2022-12-31"),
   expo = 1e5, base = "previous", rel = c(100, 90), early = 5)
ok("stkflw", adp = list(c(100, 110, 105), c(120, 118)), days = NULL, people = c(1000, 1100))
bad("stkflw", people = c(1, 2))
ok("bklog", arrive = c(0, 2, 3, 7, 9), dispose = c(5, 4, 12, NA, 15), horizon = 20, target = 2)
ok("bklog", arrive = c(1, 2, 4), dispose = c(3, NA, 6))
bad("bklog", arrive = c(1, 5), dispose = c(3, 2))
ok("trnmat", state = c("N", "N", "E", "N", "W", "W", "E", "E"), id = c(1, 1, 1, 2, 2, 2, 3, 3), time = c(1, 2, 3, 1, 2, 3, 1, 2))
ok("trnmat", state = c("a", "b", "a", "b"), id = c(1, 1, 2, 2))
bad("trnmat", state = c("a", "b"), id = 1)
ok("lftab", breaks = c(0, 3, 6, 12, 24, Inf), events = c(40, 30, 20, 10, 5), cens = c(2, 2, 3, 3, 4))
ok("lftab", breaks = c("0 to 3", "3 to 6", "6 or more"), events = c(10, 8, 4), level = 0.9, ci = "plain")
bad("lftab", breaks = 1, events = 1)

# Events in time
set.seed(36)
h1 <- sort(round(c(runif(15, 0, 100), runif(10, 40, 45)), 2))
ok("hwkfit", times = h1, hor = 100)
ok("hwkfit", times = h1, hor = 100, kernel = c("exp", "lomax"), base = "cycle", per = 50, maxit = 200)
ok("hwkfit", times = sort(c(h1, h1[3])), hor = 100)
ok("hwkfit", times = h1, hor = 100, res = 0.5, seed = 3)
bad("hwkfit", times = c(1, 2), hor = 100)
ok("trndrr", y = c(402, 377, 390, 368, 331, 310, 305), x = 2017:2023)
ok("trndrr", y = c(20, 25, 30, 33), expo = c(1000, 1500, 2500, 2600), alt = "up", level = 0.9)
ok("trndrr", y = c(5, 9, 7, 4, 6, 10, 8, 5), seas = rep(1:4, 2))
ok("trndrr", y = c(0, 0, 5))
bad("trndrr", y = c(1, 2))
ok("cusrv", y = c(10, 12, 9, 11, 10, 13, 19, 22, 25, 12), start = 6)
ok("cusrv", y = c(10, 12, 9, 11, 10, 13, 19, 22, 25, 12), expd = rep(11, 10), trans = "rossi", reset = TRUE)
bad("cusrv", y = c(1, 2))
xa <- c("2019" = 120, "2020" = 98, "2021" = 105, "2022" = 130, "2023" = 141)
xb <- c("2018" = 50, "2019" = 61, "2020" = 44, "2021" = 52, "2022" = 66, "2023" = 70)
ok("xcser", x = xa, y = xb)
ok("xcser", x = xa, y = xb, dif = TRUE)
ok("xcser", x = c(1, 3, 2, 5, 4, 6, 8, 7), y = c(2, 2, 4, 3, 6, 5, 7, 9), lagmax = 2, level = 0.9)
bad("xcser", x = c(1, 2), y = c(1, 2, 3))
bad("xcser", x = rep(1, 10), y = 1:10)
bad("xcser", x = 1:10, y = rep(0.1, 10))
ok("kmrec", time = c(3, 8, 8, 12, 20, 24), event = c(1, 1, 0, 1, 0, 1), grp = c("a", "a", "a", "b", "b", "b"))
ok("kmrec", time = c(2, 9, 30, 5, 70, 11, 12, 13), unit = rep(c("p", "q", "r"), c(3, 2, 3)), end = 100)
bad("kmrec", time = c(1, 2), unit = c("a", "b"), event = c(1, 1))
ok("hypq", lam = c(0.4, 0.7, 0.3, 0.6), pref = rbind(c(1, 2, 3), c(2, 1, 3), c(3, 2, 1), c(3, 1, 2)), mu = 1.5,
   tt = matrix(c(2, 5, 8, 6, 3, 7, 9, 6, 2, 7, 8, 3), 3))
ok("hypq", lam = c(0.2, 0.3), pref = list(c(1, 2), c(2, 1)), queue = "wait")
bad("hypq", lam = 1, pref = rbind(c(1, 1)))
ok("sdbsim", n = 6, steps = 5, a0 = 1 / 30, omega = 1 / 15, eta = 0.03, gam = 0.02, theta = 0.5, seed = 1)
ok("sdbsim", n = 5, steps = 4, mode = "pde", eps = 2, dt = 0.05, rho = matrix(0.07, 5, 5), seed = 1)
bad("sdbsim", n = 2, steps = 1, seed = 1)
bad("sdbsim", n = 6, steps = 1e9, seed = 1)
ok("sdbthr", a0 = 1, omega = 1, eps = 2, dif = 1, gam = 1, eta = 0.01)
ok("sdbthr", a0 = 1, omega = 1, eps = 0, dif = 1, gam = 1, eta = 0.01, z = 2, k = c(0.5, 1, 2))
bad("sdbthr", a0 = 0, omega = 1, eps = 0, dif = 1, gam = 1, eta = 0.01)
# a background rate and three bursts, each opened by its largest event: a series the model fits
# (events spread evenly have no aftershocks, and the fit then runs to its last step without settling)
set.seed(37)
et <- sort(round(c(runif(10, 0, 100), 20 + rexp(10, 1.5), 55 + rexp(10, 1.5), 80 + rexp(8, 1.5)), 3))
et <- et[et < 100]
em <- round(2 + (seq_along(et) %% 5) / 10, 1)
em[c(which(et >= 20)[1], which(et >= 55)[1], which(et >= 80)[1])] <- c(3.8, 4, 3.6)
# one BFGS step from a given start: the same arithmetic on every platform. A whole fit of this
# model is not: its likelihood has long flat ridges in k, c and p, and the point where the search
# stops moves with the rounding of each platform (the fits themselves are tested in R)
ok("etasf", times = et, hor = 100, start = c(0.2, 0.5, 0, 1, 1.5), maxit = 1)
ok("etasf", times = et, hor = 100, mag = em, m0 = 2, start = c(0.2, 0.5, 1, 1, 1.5), maxit = 1)
bad("etasf", times = c(1, 2), hor = 10)
mt <- sort(round(runif(30, 0, 100), 3))
ok("mvhwk", times = mt, type = rep(c("a", "b"), 15), hor = 100)
ok("mvhwk", times = mt, type = rep(c("a", "b"), 15), hor = 100, beta = 1, maxit = 100)
bad("mvhwk", times = mt, type = rep("a", 30), hor = 100, start = c(1, 2))
st <- sort(round(runif(30, 0, 100), 3))
sx <- round(runif(30, 0, 10), 3)
sy <- round(runif(30, 0, 10), 3)
ok("sthwk", times = st, x = sx, y = sy, hor = 100, win = c(0, 10, 0, 10))
ok("sthwk", times = st, x = sx, y = sy, hor = 100, par = c(0.02, 0.3, 1, 0.5), at = data.frame(t = c(50, 60), x = 5, y = 5))
bad("sthwk", times = st, x = sx, y = sy, hor = 100, win = c(0, 1, 0, 1))

# Panels: difference in differences and the event study
dd <- expand.grid(year = 2015:2020, town = LETTERS[1:8])
di <- match(dd$town, LETTERS)
dfirst <- c(2018, 2018, 2018, 2019, 2019, NA, NA, NA)
dd$bylaw <- as.integer(!is.na(dfirst[di]) & dd$year >= dfirst[di])
dd$rate <- round(30 + di + 0.5 * (dd$year - 2015) - 4 * dd$bylaw + sin(seq_len(48)), 4)
dd$first <- dfirst[di]
ok("didfe", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw)
ok("didfe", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, bacon = TRUE, loo = TRUE)
ok("didfe", y = round(dd$rate), unit = dd$town, time = dd$year, trt = dd$bylaw, fam = "pois", clus = rep(c("p", "q", "r", "s"), 12))
ok("didfe", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, mtch = 1, level = 0.9)
ok("didfe", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, mtch = 1e9, level = 0.9)
bad("didfe", y = dd$rate, unit = dd$town, time = dd$year, trt = rep(0, 48))
ok("evtst", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, leads = 2, lags = 2)
bad("evtst", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, leads = 2, lags = 2, mbar = 0, draws = 1e9,
    seed = 1)
ok("evtst", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, leads = 2, lags = 2, mbar = c(0, 1), draws = 200,
   seed = 1, bin = FALSE, ma = 2)
bad("evtst", y = dd$rate, unit = dd$town, time = dd$year, first = rep(NA, 48))
bad("evtst", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, leads = 1e9)
ok("stgdid", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, reps = 19, seed = 1)
ok("stgdid", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, meth = "dcdh", reps = 0)
ok("stgdid", y = dd$rate, unit = dd$town, time = dd$year, first = dd$first, meth = "wool", ctrl = "never", reps = 0)
bad("stgdid", y = dd$rate, unit = dd$town, time = dd$year)
sc <- expand.grid(year = 2010:2017, town = LETTERS[1:6])
si2 <- match(sc$town, LETTERS)
sc$first <- ifelse(sc$town == "A", 2015, NA)
sc$rate <- round(20 + si2 + (sc$year - 2010) * (0.2 + si2 / 20) + sin(seq_len(48)) - 3 * (sc$town == "A" & sc$year >= 2015), 4)
ok("synthc", y = sc$rate, unit = sc$town, time = sc$year, first = sc$first)
ok("synthc", y = sc$rate, unit = sc$town, time = sc$year, first = sc$first, plac = FALSE)
bad("synthc", y = sc$rate, unit = sc$town, time = sc$year, first = ifelse(si2 <= 2, 2015, NA))
ok("sdidx", y = sc$rate, unit = sc$town, time = sc$year, first = sc$first, reps = 20, seed = 1)
ok("sdidx", y = sc$rate, unit = sc$town, time = sc$year, first = ifelse(si2 <= 2, 2015, NA), level = 0.9)
bad("sdidx", y = sc$rate, unit = sc$town, time = sc$year, first = rep(NA, 48))
ok("wcbx", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw)
bad("wcbx", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, reps = 1e9, seed = 1)
ok("wcbx", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, meth = "perm", reps = 99, seed = 1)
ok("wcbx", y = dd$rate, unit = dd$town, time = dd$year, meth = "power", reps = 20, adopt = 2018, ntrt = 2, eff = c(0, 2, 5),
   seed = 1)
ok("wcbx", y = dd$rate, unit = dd$town, time = dd$year, trt = dd$bylaw, wtype = "webb", reps = 99, seed = 2)
bad("wcbx", y = dd$rate, unit = dd$town, time = dd$year)
set.seed(38)
ok("wystep", obs = c(2.5, 0.4, -2.1), null = matrix(round(rnorm(90), 3), ncol = 3))
bad("wystep", obs = c(1, 2), null = matrix(1, 2, 3))

# Covariate balance: effects
set.seed(39)
cx <- matrix(round(rnorm(80), 3), 40)
cd <- rbinom(40, 1, plogis(cx[, 1]))
cy <- round(2 * cd + cx[, 1] - cx[, 2] + rnorm(40), 3)
ok("cbatt", y = cy, d = cd, x = cx, reps = 19, seed = 1)
ok("cbatt", y = cy, d = cd, x = data.frame(a = cx[, 1], b = cx[, 2]), est = "ebal", tgt = "att", reps = 0)
ok("cbatt", y = cy, d = cd, x = cx, est = "psm", tgt = "att", calip = 0.5)
ok("cbatt", y = cy, d = cd, x = cx, est = "gcomp", reps = 0, level = 0.9)
bad("cbatt", y = cy, d = cd, x = cx, est = "ebal", tgt = "ate")

# Surveys
set.seed(40)
sv <- data.frame(region = rep(c("north", "south"), each = 20), block = rep(1:8, each = 5),
                 victim = rbinom(40, 1, rep(c(0.2, 0.35), each = 20)), w = round(runif(40, 50, 150)))
ok("svymx", y = sv$victim, w = sv$w, strata = sv$region, psu = sv$block, stat = "prop")
ok("svymx", y = sv$victim, w = sv$w, strata = sv$region, psu = sv$block, stat = "prop", by = sv$region, ref = "south")
ok("svymx", y = sv$w, den = rep(2, 40), stat = "ratio", fpc = rep(1000, 40), level = 0.9)
bad("svymx", y = c(1, 2), stat = "prop")

# Learners: double machine learning and the super learner
set.seed(41)
lx <- matrix(round(rnorm(60 * 3), 3), 60)
ld <- round(lx[, 1] + rnorm(60), 3)
ly <- round(0.5 * ld + lx[, 1] - lx[, 2] + rnorm(60), 3)
lf <- sample(rep(1:3, length.out = 60))
ok("dmlx", y = ly, d = ld, x = lx, fold = lf, lrn = "glm", lrnd = "glm")
bad("dmlx", y = ly, d = ld, x = lx, fold = lf, lrn = "glm", lrnd = "glm", nlam = 1e9)
ok("dmlx", y = ly, d = ld, x = lx, fold = lf, lrn = "glm", lrnd = "lasso", score = "iv", nlam = 10)
ok("dmlx", y = ly, d = ld, x = lx, fold = cbind(lf, rev(lf)), lrn = "lasso", lrnd = "glm", nlam = 10)
bad("dmlx", y = ly, d = ld, x = lx, fold = rep(1, 60))
ok("suplrn", x = lx, y = ly, lrns = c("glm", "lasso"), nfold = 3, nlam = 10)
bad("suplrn", x = lx, y = ly, lrns = c("glm", "lasso"), nfold = 3, nlam = 1e9)
ok("suplrn", x = lx, y = as.integer(ly > 0), fam = "binomial", lrns = c("glm", "forest"), fold = lf, ntree = 20, newx = lx[1:3, ],
   seed = 2)
bad("suplrn", x = lx, y = ly, lrns = c("glm", "glm"))

# Spatial methods
set.seed(42)
sp <- cbind(x = round(runif(40, 0, 2000)), y = round(runif(40, 0, 2000)))
ok("knxnr", x = sp, t = round(runif(40, 0, 365)), ds = c(200, 400), dt = c(14, 28), nsim = 19, seed = 2)
ok("knxnr", x = sp, t = round(runif(40, 0, 365)), ds = 300, dt = 30, nsim = 0)
bad("knxnr", x = sp[1:2, ], t = c(1, 2), ds = 1, dt = 1)
ok("knxnr", x = sp, t = 1:40, ds = 5e-324, dt = 30, nsim = 0)
kp <- cbind(x = round(c(rnorm(15, 2, 0.5), runif(10, 0, 6)), 3), y = round(c(rnorm(15, 3, 0.5), runif(10, 0, 6)), 3))
ok("kdsurf", x = kp, win = c(0, 6, 0, 6), nx = 5, edge = "uniform")
ok("kdsurf", x = kp, h = 0.4, at = "points")
ok("kdsurf", x = kp, h = c(0.5, 0.7), at = cbind(c(1, 2), c(3, 4)), w = rep(c(1, 2), length.out = 25))
bad("kdsurf", x = kp, edge = "diggle")
dp <- rbind(cbind(x = round(rnorm(10, 0, 0.2), 3), y = round(rnorm(10, 0, 0.2), 3)),
            cbind(x = round(rnorm(10, 3, 0.2), 3), y = round(rnorm(10, 1, 0.2), 3)), cbind(x = c(-2, 5), y = c(3, -2)))
ok("dbscn", x = dp, eps = 0.3, minpts = 4)
ok("dbscn", x = dp, eps = 0.5, minpts = 3, border = FALSE)
bad("dbscn", x = dp, eps = 0)
ok("dbscn", x = dp, eps = 5e-324, minpts = 2)
hp <- cbind(x = round(c(rnorm(20, 2, 0.3), runif(20, 0, 10)), 3), y = round(c(rnorm(20, 2, 0.3), runif(20, 0, 10)), 3))
ok("hsper", x = hp, prd = rep(2021:2024, length.out = 40), cell = 1, top = 0.1)
ok("hsper", x = hp, prd = rep(1:2, 20), cell = 2, org = c(-1, -1))
bad("hsper", x = hp, prd = rep(1:2, 20), cell = 1, top = 2)
ok("ebshr", y = c(40, 12, 9, 25, 7, 31, 5, 18))
ok("ebshr", y = c(40, 12, 9, 25, 7, 31, 5, 18), v = 10, w = c(1, 2, 1, 2, 1, 2, 1, 2), b = 0.5)
bad("ebshr", y = c(-5, -3))
kt <- c(sample(0:99, 30, TRUE), sample(40:44, 10, TRUE))
kq <- cbind(x = round(c(runif(30, 0, 20), rnorm(10, 5, 0.3)), 3), y = round(c(runif(30, 0, 20), rnorm(10, 15, 0.3)), 3))
ok("kscan", x = kq, t = kt, r = c(1, 3), tw = c(5, 20), top = 2, nsim = 19, seed = 1)
ok("kscan", x = kq, t = kt, r = 2, tw = 10, nsim = 0, minc = 3)
bad("kscan", x = kq[1:2, ], t = c(1, 2), r = 1, tw = 1)
gp <- cbind(x = round(runif(20), 3), y = round(runif(20), 3))
gw <- spwt(gp, k = 3)
ok("gmoran", x = round(gp[, 1] + rnorm(20, sd = 0.2), 3), w = gw, nsim = 19, seed = 1)
ok("gmoran", x = data.frame(y1 = rpois(20, 3), y2 = rpois(20, 3 + 4 * gp[, 1])), w = gw)
bad("gmoran", x = c(1, 2, 3), w = gw)
ok("bvmor", x = round(gp[, 1] + rnorm(20, sd = 0.3), 3), y = round(gp[, 1] + rnorm(20, sd = 0.3), 3), w = gw, nsim = 19, seed = 1)
ok("bvmor", x = round(gp[, 2], 3), y = round(gp[, 1], 3), w = gw, nsim = 0)
bad("bvmor", x = rep(1, 20), y = gp[, 1], w = gw, nsim = 0)
sW <- matrix(0, 20, 20)
sW[cbind(gw$i, gw$j)] <- gw$w
sxv <- round(rnorm(20), 3)
syv <- round(as.vector(solve(diag(20) - 0.4 * sW, 1 + 2 * sxv + rnorm(20))), 3)
ok("sarml", y = syv, x = sxv, w = gw)
ok("sarml", y = syv, x = data.frame(a = sxv, b = round(gp[, 1], 3)), w = gw, model = "sac", durbin = TRUE)
ok("sarml", y = syv, x = sxv, w = sW, method = "2sls")
bad("sarml", y = syv[1:3], x = sxv[1:3], w = gw)
ok("spill", w = gw, rho = 0.4, b = c(gent = 0.1), th = 0.3)
ok("spill", fit = sarml(syv, sxv, gw))
bad("spill", w = gw, rho = 0.4)
ck <- cbind(x = round(runif(30), 3), y = round(runif(30), 3))
cm <- ifelse(ck[, 1] < 0.5, "burglary", "robbery")
ok("crossk", x = ck, marks = cm, i = "burglary", j = "robbery", win = c(0, 1, 0, 1), r = c(0.05, 0.1, 0.2))
ok("crossk", x = ck, marks = cm, i = 1, j = 2, win = c(0, 1, 0, 1), r = c(0.1, 0.2), sigma = 0.1, nsim = 9, seed = 1)
bad("crossk", x = ck, marks = cm, i = "theft", j = "robbery", win = c(0, 1, 0, 1), r = 0.1)
qp <- cbind(x = round(runif(30, 0, 2), 3), y = round(runif(30), 3))
ok("qdrtst", x = qp, win = c(0, 2, 0, 1), nx = 4, ny = 2)
ok("qdrtst", x = qp, win = c(0, 2, 0, 1), nx = 3)
bad("qdrtst", x = qp, win = c(0, 2, 0, 1), nx = 1, ny = 1)
ok("morsta", counts = c(3, 0, 7, 12, 1, 4, 0, 9))
ok("morsta", counts = rbind(c(3, 0), c(7, 12)), crit = 0.1)
bad("morsta", counts = c(1, 0))
dx <- cbind(a = c(90, 80, 10, 5, 50), b = c(10, 20, 85, 95, 50))
ok("dissix", x = unname(dx))
ok("dissix", x = data.frame(dx), delta = 0.2, w = spwt(nb = list(2, c(1, 3), c(2, 4), c(3, 5), 4), type = "contig"))
bad("dissix", x = cbind(c(1, 2)))
sq <- function(x0, id) data.frame(x = x0 + c(0, 1, 1, 0), y = c(0, 0, 1, 1), id = id)
pp2 <- rbind(sq(0, "A"), sq(1, "B"), sq(2, "C"))
ok("pinpl", x = cbind(c(0.5, 1.5, 2.5, 4), c(0.5, 0.5, 0.2, 0.5)), poly = pp2)
ok("pinpl", x = cbind(c(0.5, 1.5, 0.2), c(0.5, 0.5, 0.9)), poly = pp2, pop = c(100, 200, 300))
bad("pinpl", x = cbind(0.5, 0.5), poly = list(a = 1))
ok("arwt", x = c(a = 120, b = 40, c = 60, d = 100), from = c("a", "b", "c", "d", "d"), to = c("N1", "N1", "N2", "N2", "N3"),
   area = c(1, 1, 1, 7, 3), top = 1)
ok("arwt", x = c(5, 6), id = c("p", "q"), from = c("p", "q"), to = c("Z", "Z"))
bad("arwt", x = c(1, 2), from = "a", to = "b")
set.seed(43)
bp <- cbind(x = round(runif(30, 0, 5000)), y = round(runif(30, 0, 5000)))
ok("bufas", x = bp, site = cbind(c(1000, 3500), c(1200, 3000)), comp = cbind(c(2500, 4500), c(500, 4000)), r = 600)
ok("bufas", x = bp, site = cbind(1000, 1200), comp = cbind(2500, 500), r = 900, buf = 1000)
bad("bufas", x = bp, site = cbind(1000, 1200), comp = cbind(2500, 500), r = 0)

# clprep(): one case for every preset, on a table of eight rows made from
# the preset's own rows (the labels its recodes name, numbers where it
# reads numbers, one missing value), then presets of one's own for the
# recipes and recodes the shipped presets do not use, and the refusals
pt <- clprep()
tb <- utils::read.delim(system.file("presets", "presets.tsv", package = "crimlab"),
                        colClasses = "character", na.strings = character(0), quote = "")
made <- list(long = c("count", "with"), acvm = c("ac", "vm", "hiac"), pyr = c("np", "nreg", "ac"),
             acnv = c("acn", "vmany", "hiacn"), pairb = c("tb", "yb"), pairc = c("tc", "yc"),
             gent = "gent", pop = "pop")
otc <- c("UniqueIndividual_ID", "EndFiscalYear", "Region_AtTimeOfPlacement",
         "Region_MostRecentPlacement", "MentalHealth_Alert", "SuicideRisk_Alert", "SuicideWatch_Alert")
cap1 <- function(s) paste0(toupper(substr(s, 1L, 1L)), substring(s, 2L))
fixed <- list(
  UniqueIndividual_ID = c(101, 101, 102, 101, 103, 102, 101, 103),
  EndFiscalYear = c(2020, 2020, 2020, 2021, 2021, 2020, 2020, 2021),
  Region_AtTimeOfPlacement = c("East", "East", "North", "West", "East", "North", "Central", "West"),
  Region_MostRecentPlacement = c("East", "West", "North", "West", "Central", "North", "Central", "East"),
  MentalHealth_Alert = c("Yes", "No", "yes", "No", "Y", "no", "TRUE", "No"),
  SuicideRisk_Alert = c("No", "Yes", "no", "No", "Yes", "no", "1", "No"),
  SuicideWatch_Alert = c("No", "No", "Yes", "No", "yes", "No", "No", "Yes"),
  Gender = c("Female", "Male", "Male", "Female", "Male", "Female", "Male", "Male"),
  Age_Category = c("18 to 49", "50 and over", "18 to 49", "18 to 49", "50 and over", "18 to 49",
                   "50 and over", "18 to 49"),
  Number_Of_Placements = c(1, 2, 3, 1, 2, 1, 4, 1),
  Measure = c("Median", "Median", "Max", "Median", "Max", "Median", "Median", "Max")
)
num8 <- c(12, 30, 8, 25, 14, 41, 9, 22)
synv <- function(col, rec) {
  if (!is.null(fixed[[col]])) return(fixed[[col]])
  map <- rec[!rec %in% c("", "fac", "num", "pct", "one", "flag", "dum") & !startsWith(rec, "sub:") &
               !startsWith(rec, "const:")]
  if (length(map)) {
    rules <- unlist(strsplit(map, ";", fixed = TRUE))
    keys <- vapply(strsplit(rules, "=", fixed = TRUE), `[`, "", 1L)
    other <- all(grepl("(^|;)\\*=", map))
    keys <- keys[keys != "*"]
    labs <- unlist(strsplit(keys, "|", fixed = TRUE))
    labs <- ifelse(startsWith(labs, "has:"), paste("with", substring(labs, 5L)), labs)
    labs <- cap1(unique(c(labs, if (other) "other")))
    return(rep(labs, length.out = 8L))
  }
  if (any(startsWith(rec, "sub:"))) {
    return(rep(c("2019-04-01", "2020-03-15", "2021-11-30"), length.out = 8L))
  }
  if (any(rec == "flag")) return(c("Yes", "no", "TRUE", "0", "y", "No", "t", "no"))
  if (any(rec %in% c("fac", "dum")) && !any(rec %in% c("", "num", "pct"))) {
    return(rep(c("North", "East", "West"), length.out = 8L))
  }
  num8
}
syn <- function(r) {
  rec <- r$arg == "recipe"
  out <- unique(unlist(made[r$column[rec]]))
  need <- character(0)
  for (i in which(rec)) {
    a <- r$recode[i]
    need <- c(need, switch(r$column[i],
      long = strsplit(a, ",", fixed = TRUE)[[1L]],
      keep = strsplit(a, ",", fixed = TRUE)[[1L]],
      acvm = , acnv = c(otc, "Gender", "Age_Category"),
      pyr = otc,
      pairb = c(otc[-(3:4)], "Number_Of_Placements"),
      pairc = c(otc[3:4], "NumberConsecutiveDays_Segregation"),
      pop = strsplit(a, ";", fixed = TRUE)[[1L]][1:2],
      gent = strsplit(a, ",", fixed = TRUE)[[1L]][1:4]))
  }
  fc <- as.character(unlist(strsplit(r$filter[nzchar(r$filter)], ";", fixed = TRUE)))
  fc <- vapply(strsplit(fc, "=", fixed = TRUE), `[`, "", 1L)
  cols <- unique(c(need, setdiff(c(r$column[!rec & nzchar(r$column)], fc), out)))
  d <- as.data.frame(lapply(stats::setNames(cols, cols), function(k) {
    synv(k, r$recode[!rec & r$column == k])
  }), check.names = FALSE)
  d[3L, cols[1L]] <- NA
  d
}
for (p in pt$preset) ok("clprep", data = syn(tb[tb$preset == p, ]), preset = p)
ok("clprep")
g6 <- data.frame(inc = c(10, 30, 12, 50, 9, 20), rent = c(5, 9, 4, 12, 6, 5),
                 colg = c(0.3, 0.1, 0.5, 0.2, 0.6, 0.1), rentg = c(0.2, 0.4, 0.5, 0.1, 0.3, NA),
                 y = c(3, 1, 4, 1, 5, 9))
ok("clprep", data = g6, preset = data.frame(
  preset = "mine", generic = "cntreg", arg = c("recipe", "y", "grp"),
  column = c("gent", "y", "gent"), recode = c("inc,rent,colg,rentg", "", "fac")))
ok("clprep", data = g6, preset = data.frame(
  preset = "mine", generic = "cntreg", arg = c("recipe", "y", "grp"),
  column = c("gent", "y", "gent"), recode = c("inc,rent,colg,rentg,0.4,0.5", "", "")))
ok("clprep", data = data.frame(a = c("x", "1", "2.5"), b = c(4, 5, 6), w = c("2020-01-02", "2021-05-06", NA)),
   preset = data.frame(preset = "mine", generic = "grprr", arg = c("x", "pop", "lv", "k", "yr", "yr"),
                       column = c("a", "b", "", "", "w", "b"),
                       recode = c("num", "pct", "const:TRUE", "const:per", "sub:3,4", "const:2.5"),
                       units = c("", "people", "", "", "", ""), note = c("first", "", "", "", "", "second")))
ok("clprep", data = data.frame(a = c(1, 2), b = c(3, 4), k = c("p", "q")),
   preset = data.frame(preset = "mine", generic = "cntreg", arg = c("recipe", "y", "grp", "adj"),
                       column = c("long", "count", "with", "k"), recode = c("a,b", "", "", "dum")))
ok("clprep", data = data.frame(n = c(100, 30), r = c("50", "x"), z = c(" A ", "b")),
   preset = data.frame(preset = "mine", generic = "grprr", arg = c("recipe", "x", "g", "h"),
                       column = c("pop", "pop", "z", "z"),
                       recode = c("n;r;1000", "", "a=1;has:b=2", "has:=7"), filter = c("", "", "", "")))
bad("clprep", data = data.frame(a = 1))
bad("clprep", data = data.frame(a = 1), preset = c("a", "b"))
bad("clprep", data = data.frame(a = 1), preset = "no-such-preset")
bad("clprep", data = list(1, 2), preset = "otis-c01-cnt")
bad("clprep", data = data.frame(a = 1), preset = "otis-c01-cnt")
bad("clprep", data = data.frame(s = c("x", "y")),
    preset = data.frame(preset = "m", generic = "g", arg = "a", column = "s", recode = "x=1"))
bad("clprep", data = data.frame(s = 1),
    preset = data.frame(preset = "m", generic = "g", arg = "recipe", column = "melt"))
bad("clprep", data = data.frame(s = 1), preset = data.frame(preset = "m", generic = "g"))
bad("clprep", data = data.frame(s = 1), preset = data.frame(preset = c("m", "n"), generic = "g", arg = "a"))
bad("clprep", data = data.frame(s = 1),
    preset = data.frame(preset = "m", generic = "g", arg = "recipe", column = "keep", recode = "s,t"))

# Pedigree likelihood ratios
tri <- data.frame(id = c("AF", "M", "C"), fid = c(NA, NA, "AF"), mid = c(NA, NA, "M"))
tra <- data.frame(id = c("AF", "M", "C", "X"), fid = c(NA, NA, "X", NA), mid = c(NA, NA, "M", NA))
tg <- data.frame(id = rep(c("AF", "M", "C"), 2), locus = rep(c("D3", "vWA"), each = 3),
                 a1 = c(14, 15, 15, 16, 17, 17), a2 = c(16, 17, 16, 18, 17, 19))
tf <- data.frame(locus = rep(c("D3", "vWA"), c(4, 4)), allele = c(14, 15, 16, 17, 16, 17, 18, 19),
                 freq = c(0.1, 0.25, 0.25, 0.2, 0.2, 0.3, 0.2, 0.1))
ok("pedlr", ped = tri, geno = tg, freq = tf, alt = tra)
ok("pedlr", ped = tri, geno = tg, freq = tf, alt = tra, mut = "equal", rate = c(0.001, 0.003))
sb <- data.frame(id = c("f", "m", "a", "b"), fid = c(0, 0, "f", "f"), mid = c(0, 0, "m", "m"))
ok("pedlr", ped = sb, geno = data.frame(id = c("a", "b"), locus = "L", a1 = c(9, 11), a2 = c(12, 13)),
   freq = data.frame(locus = "L", allele = 9:13, freq = c(0.1, 0.2, 0.3, 0.25, 0.15)),
   mut = "step", rate = 0.005, range = 0.4)
bad("pedlr", ped = tri, geno = tg, freq = tf, mut = "stepwise")
bad("pedlr", ped = tri[2:3, ], geno = tg, freq = tf)
bad("pedlr", ped = tri, geno = tg, freq = tf, mut = "step")
# The continuous mixture model of peak heights
gp <- data.frame(locus = rep(c("D3", "vWA"), c(4, 3)), allele = c(14, 15, 16, 17, 16, 17, 18),
                 height = c(1210, 380, 1150, 420, 1400, 1960, 610))
gf <- data.frame(locus = rep(c("D3", "vWA"), c(5, 5)), allele = c(13, 14, 15, 16, 17, 15, 16, 17, 18, 19),
                 freq = c(0.05, 0.1, 0.25, 0.25, 0.2, 0.1, 0.2, 0.3, 0.2, 0.1),
                 size = c(110, 114, 118, 122, 126, 160, 164, 168, 172, 176))
gr <- data.frame(id = rep(c("V", "S"), each = 2), locus = rep(c("D3", "vWA"), 2),
                 a1 = c(14, 16, 15, 17), a2 = c(16, 17, 17, 18))
ok("gmixlr", peaks = gp, freq = gf, ref = gr, hp = c("V", "S"), nc = 2, hd = "V",
   par = c(mu = 1800, omega = 0.1, m1 = 0.73, m2 = 0.27))
ok("gmixlr", peaks = gp, freq = gf, ref = gr, hp = c("V", "S"), nc = 2, hd = "V", theta = 0.03,
   stut = TRUE, deg = TRUE, par = c(mu = 1900, omega = 0.12, beta = 0.9, xi = 0.05, m1 = 0.7, m2 = 0.3))
ok("gmixlr", peaks = gp[c(1, 3, 4, 5, 6), ], freq = gf, ref = gr, hp = "S", nc = 2, at = 100,
   pc = 0.1, lambda = 0.02, par = c(mu = 1500, omega = 0.2, m1 = 0.25, m2 = 0.75))
# no case without `par`: a free fit ends after a number of steps that
# differs between machines (9 on x86-64, 13 on arm64 for these data), and
# `iter` is part of the result; the fit is tested in the R package against
# EuroForMix (tests/testthat/test-gmixlr.R)
bad("gmixlr", peaks = gp, freq = gf, ref = gr, hp = "S", nc = 0)
bad("gmixlr", peaks = gp, freq = gf, ref = gr, hp = "W", nc = 2)
bad("gmixlr", peaks = gp, freq = gf[1:3], ref = gr, hp = "S", nc = 1, deg = TRUE)
# The CMPS score of bullet signatures
ca <- round(cumsum(sin(seq_len(400) / 3) + cos(seq_len(400) / 7)) + sin(seq_len(400) * 1.3), 4)
cb <- c(ca[26:400], rep(NA, 5)) + round(0.05 * cos(seq_len(380) * 2.1), 4)
ok("cmpsc", a = c(NA, NA, ca[-(1:2)]), b = cb)
ok("cmpsc", a = ca, b = cb, len = 40, npeaks = c(4, 2), tx = 10)
ok("cmpsc", a = ca, b = rev(cb), npeaks = 3)
ok("cmpsc", a = ca, b = cb, len = 40, full = TRUE)
bad("cmpsc", a = ca, b = cb, npeaks = c(1, 0))
bad("cmpsc", a = ca[1:20], b = cb)


# Discrimination in searches
tsd <- data.frame(
  grp = rep(c("a", "b"), 4), area = rep(1:4, each = 2),
  stops = c(1200, 800, 900, 1100, 1500, 700, 1000, 1000),
  srch = c(60, 70, 50, 95, 80, 66, 52, 90),
  hits = c(24, 18, 19, 25, 30, 17, 21, 24)
)
ok("thrtst", stops = c(100, 100, 80), srch = c(20, 40, 16), hits = c(8, 10, 4), grp = c("a", "b", "c"), ref = "b")
ok("thrtst", stops = tsd$stops, srch = tsd$srch, hits = tsd$hits, grp = tsd$grp, area = tsd$area, maxit = 1)
ok("thrtst", stops = c(100, 100, 100, 100), srch = c(10, 20, 30, 10), hits = c(5, 5, 5, 5), grp = c(1, 1, 2, 2),
   off = c("x", "y", "x", "y"), level = 0.9)
ok("thrtst", stops = rep(600, 9), srch = c(72, 48, 60, 60, 72, 48, 48, 60, 72), hits = rep(c(36, 24, 30), 3),
   grp = rep(1:3, each = 3), off = rep(c("a", "b", "c"), 3), rank = "joint")
bad("thrtst", stops = c(10, 10), srch = c(11, 1), hits = c(1, 1), grp = c(1, 2))
tbd <- expand.grid(g = c("a", "b", "c"), a = 1:5, stringsAsFactors = FALSE)
tbs <- c(60, 50, 40, 55, 45)[tbd$a]
tbk <- c(12, 9, 8, 11, 7)[tbd$a]
tbh <- c(5, 3, 3, 4, 2)[tbd$a]

# Internal benchmarking of officers
set.seed(5)
ibo <- rep(c("a", "b", "c"), each = 10)
ibx <- rbinom(30, 1, ifelse(ibo == "a", 0.7, 0.3))
iby <- c(1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0)
ok("ibench", y = iby, off = ibo, x = data.frame(night = ibx), minn = 5)
ok("ibench", y = iby, off = ibo, x = data.frame(night = ibx), minn = 5, dr = TRUE)
bad("ibench", y = iby, off = ibo, minn = 5, local = TRUE)
ok("ibench", y = c(1, 0, 1, 1, 0, 0, 1, 0, 1), off = c("a", "a", "a", "a", "b", "b", "b", "b", "c"), minn = 2,
   fdr = 0.1)
bad("ibench", y = c(1, 0), off = c("a", "a"))

# Eyewitness lineups
set.seed(3)
ltp <- rep(0:1, each = 20)
lpr <- rep(c("sim", "seq"), 20)
lpk <- ifelse(ltp == 1, sample(c("suspect", "filler", "none"), 40, TRUE, c(0.5, 0.15, 0.35)),
              sample(c("suspect", "filler", "none"), 40, TRUE, c(0.15, 0.3, 0.55)))
lcf <- ifelse(lpk == "none", NA, sample(1:5, 40, TRUE))
ok("lnroc", tp = ltp, pick = lpk, conf = lcf, proc = lpr)
ok("lnroc", tp = ltp, pick = lpk, conf = lcf, proc = lpr, far = 0.1, nboot = 20, seed = 1, level = 0.9)
ok("lnroc", tp = c(1, 1, 1, 1, 0, 0, 0, 0), pick = c("suspect", "suspect", "filler", "none", "suspect", "filler",
   "none", "none"), conf = c(3, 1, 2, NA, 2, 3, NA, NA), size = 3)
bad("lnroc", tp = c(1, 0), pick = c("suspect", "guess"), conf = c(1, 1))
ok("lnsize", choice = c(rep(3, 15), rep(1, 5), rep(2, 6), rep(4, 4), rep(5, 7), rep(6, 3)), k = 6, susp = 3)
ok("lnsize", choice = c(1, 1, 1, 2, 3, 3), level = 0.9)
bad("lnsize", choice = c(1, 1, 1, 2, 3, 3), k = 1e9)
ok("lnsize", choice = c(0, 2, 2, 2, 3, 3, 4), k = 4, susp = 2)
ok("lnroc", tp = ltp, pick = lpk, conf = lcf, proc = lpr, far = 0.5)
bad("lnsize", choice = c(1, 7), k = 6)

# Repeat victimisation
ok("rptvic", id = c(1, 2, 1, 3, 1, 2), t = c(1, 5, 4, 9, 20, 6), pop = 10, brk = c(0, 7, 30))
ok("rptvic", id = c("a", "b", "a", "c", "c", "c", "d"))
ok("rptvic", id = c("a", "b", "a", "b"), t = c(3, 1, 10, 2))
bad("rptvic", id = c("a", "b"), pop = 1)

# The weighted displacement quotient
ok("wdqx", a = c(100, 60), b = c(50, 40), c = c(200, 180))
ok("wdqx", a = c(10, 10), b = c(5, 9), c = c(10, 10))
ok("wdqx", a = c(207, 110), b = c(178, 157), c = c(308, 318))
bad("wdqx", a = c(1, 2), b = c(1, 2), c = c(0, 2))

# Prospective hot-spotting
set.seed(6)
php <- cbind(x = round(runif(40, 0, 500)), y = round(runif(40, 0, 400)))
phd <- sample(1:60, 40, TRUE)
ok("prohs", x = php, t = phd, at = 45, win = c(0, 500, 0, 400), cell = 100, bw = 250, tw = 30, ahead = 14)
ok("prohs", x = cbind(c(0.5, 0.6, 1.6), c(0.5, 0.4, 0.5)), t = c(0, 8, 9), at = 7, win = c(0, 2, 0, 1), cell = 1,
   bw = 2, top = 0.5)
bad("prohs", x = cbind(c(1, 2), c(1, 2)), t = c(1, 2), at = 2, win = c(0, 5, 0, 5), cell = 1, top = 0)

# Crime linkage
lkx <- cbind(c(20, 20, 49, 7, -8, 26, -1, 4, 25, 55, 14, 29), c(45, 60, 41, 41, 21, 49, 28, 45, 46, 35, 22, 45))
lkt <- c(35, 59, 27, 70, 30, 26, 54, 30, 5, 18, 8, -6)
ok("lnkpr", x = lkx, t = lkt, series = c(1, 1, 1, 1, NA, 2, 2, 2, 3, 3, NA, 3))
ok("lnkpr", x = lkx, t = lkt, series = c("a", "a", "a", "a", NA, "b", "b", "b", "c", "c", NA, "c"),
   mo = cbind(c(1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1), c(0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0)))
bad("lnkpr", x = lkx[1:3, ], t = 1:3, series = c(1, 1, 1))

# The split-population model, from its fixed start (a few steps: the
# stopping point of a converged fit moves with rounding)
set.seed(10)
scx <- rbinom(30, 1, 0.4)
sct <- round(rweibull(30, 1.3, 30), 2)
sce <- as.numeric(sct < 40 & runif(30) < 0.7)
sct[sce == 0] <- pmax(sct[sce == 0], 40)
ok("spcure", time = sct, event = sce, x = data.frame(prior = scx), maxit = 3)
ok("spcure", time = sct, event = sce, dist = "lognormal", maxit = 3, level = 0.9)
ok("spcure", time = c(1, 2, 3, 4, 50, 60, 70, 80), event = c(1, 1, 1, 1, 0, 0, 0, 0), dist = "exp", maxit = 2)
bad("spcure", time = c(1, 2), event = c(1, 1))

# The age-period-cohort model
set.seed(11)
acy <- matrix(rpois(20, 50 * exp(seq(-0.5, 0.5, length.out = 20))), 5)
acp <- matrix(round(runif(20, 800, 1200)), 5)
ok("apcie", y = acy, expo = acp)
ok("apcie", y = matrix(c(5, 8, 3, 9, 12, 6, 4, 10, 7), 3), expo = matrix(100, 3, 3), level = 0.9)
bad("apcie", y = matrix(1, 2, 3), expo = matrix(1, 2, 3))

# The accuracy of a hot-spot map
ok("hsacc", pred = c(9, 8, 7, 6, 5, 4, 3, 2, 1, 0), obs = c(3, 0, 5, 1, 0, 0, 2, 9, 0, 0), top = c(0.2, 0.3),
   prev = c(4, 2, 2, 0, 0, 0, 1, 1, 0, 0))
ok("hsacc", pred = c(FALSE, TRUE, FALSE, FALSE), obs = c(1, 4, 0, 2), area = c(1, 2, 1, 1))
bad("hsacc", pred = c(1, 2), obs = c(0, 0))
ok("hsacc", pred = c(5, 5, 5, 5, 4, 3, 2, 1, 0, 0), obs = c(3, 0, 5, 1, 0, 0, 2, 9, 0, 0), top = 0.2, ties = "frac")
ok("hsacc", pred = c(5, 5, 5, 5, 4, 3, 2, 1, 0, 0), obs = c(3, 0, 5, 1, 0, 0, 2, 9, 0, 0), top = 0.2, ties = "out")
bad("hsacc", pred = c(1, 2), obs = c(1, 0), ties = "all")

# Bounds under monotone response and selection
ok("mtrbnd", y = c(1, 0, 0, 1, 1, 0, 1, 1), z = c(1, 1, 1, 2, 2, 2, 3, 3))
ok("mtrbnd", y = c(2, 5, 3, 9), z = c("a", "a", "b", "b"), ylo = 0, yhi = 10, w = c(1, 3, 1, 1), assume = "mtr")
ok("mtrbnd", y = c(1, 0, 0, 1, 1, 0), z = c(1, 1, 2, 2, 3, 3), assume = "mts")
bad("mtrbnd", y = c(0, 2), z = c(1, 2))

# Search theory for search and rescue
ok("ltrlrng", x = c(-3, -0.7, 0, 0.2, 0.5, 1, 2.5), W = 1.3)
ok("ltrlrng", x = c(-1, 0, 0.3, 2), W = 1, model = "mbeta", m = 0.4, navsd = 0.25)
ok("ltrlrng", x = c(0, 0.5, 1), W = 0.8, navsd = 0.3)
bad("ltrlrng", x = 1, W = -1)
bad("ltrlrng", x = 1, W = 1, model = "bogus")
ok("swpwdth", x = seq(-2, 2, by = 0.5), p = c(0, 0.1, 0.3, 0.8, 1, 0.8, 0.3, 0.1, 0))
bad("swpwdth", x = c(2, 1), p = c(0, 1))
ok("srchpod", coverage = c(0, 0.5, 1, 2), model = "invcube")
ok("srchpod", coverage = c(0.2, 0.9, 1.7), model = "mbeta", m = 0.5)
ok("srchpod", W = 0.5, effort = c(1, 2, 3), area = c(2, 2, 4))
bad("srchpod", coverage = 1, W = 0.6)
ok("pswppod", W = 1, S = 1.3, npts = 101, ntrk = 20)
ok("pswppod", W = 1, S = 0.6, model = "mbeta", m = 0.8, navsd = 0.25, npts = 51, ntrk = 10)
bad("pswppod", W = 1, S = 1, npts = 3.5)
ok("srchscs", poc = c(0.4, 0.3, 0.2, 0.05), pod = rbind(c(0.8, 0.5, 0, 0.3), c(0.2, 0.6, 0.9, 0.1)))
ok("srchscs", poc = 1, pod = 1)
cm <- datumpmap(seq(-2, 2), seq(-2, 2), sigma = 1.5)
ok("srchscs", poc = cm, pod = matrix(0.5, 4, 4))
bad("srchscs", poc = c(0.8, 0.8), pod = c(0.5, 0.5))
ok("datumpmap", xedg = seq(-2, 2), yedg = seq(-1, 2), sigma = 1.5, center = c(0.3, 0.2))
ok("datumpmap", xedg = c(-2, -0.5, 1, 3), yedg = c(-1, 0, 2), datum = "line", center = c(-1, -0.5),
   end = c(2, 1), nodes = 16)
ok("datumpmap", xedg = 0:4, yedg = 0:2, datum = "area", center = c(1, 0.5), end = c(3, 1.5))
ok("datumpmap", xedg = 0:4, yedg = c(-1.5, -0.5, 0.5, 1.5), datum = "corridor", center = c(0.5, 0), end = c(3, 0))
ok("datumpmap", xedg = 0:3, yedg = 0:2, datum = "general", wts = matrix(1:6, 2))
bad("datumpmap", xedg = 1:2, yedg = 1:2, datum = "line")
ok("optcrcsrch", sigma = 2, W = 1.5, effort = 30)
bad("optcrcsrch", sigma = 1, W = 1, effort = -1)
ok("opteffort", poc = c(0.5, 0.2, 0.1), area = c(1, 2, 0.5), W = 1, effort = 3)
ok("opteffort", poc = cm, W = 0.8, effort = 6)
bad("opteffort", poc = c(0.5, 0.5), W = 1, effort = 1)
ok("optsqr", effort = c(4, 4, 8), E = 2)
ok("optsqr", effort = 4, model = "mbeta", m = 0.5)
bad("optsqr", effort = -1)
ok("srchfig")
ok("srchfig", fig = 5, num = 8)
bad("srchfig", fig = "9-9")

# Air pollution
arc <- rbind(c(300, 20, 0), c(1500, -40, 2), c(5000, 0, 10), c(-1, 0, 0))
ok("pgsig", x = c(50, 400, 3000), stab = "F", set = "urban")
bad("pgsig", x = 100, stab = "G")
ok("briggs", x = c(10, 500, 1e5), u = 4, diam = 1, exitvel = 12, stktmp = 420, ambtmp = 285)
ok("briggs", x = c(10, 1e5), u = 4, diam = 1, exitvel = 12, stktmp = 420, ambtmp = 285, stab = "F", dthdz = 0.01)
bad("briggs", x = 10, u = 4, diam = 1, exitvel = 12, stktmp = 155, ambtmp = 20)
ok("gplume", q = 10, u = 3, h = 50, rcpt = arc, stab = "A", set = "urban", mixh = 300, nimg = 2)
ok("gplume", q = 10, u = 3, h = 0, rcpt = c(800, 0, -1, 900, 5, 0))
bad("gplume", q = 10, u = 3, h = 0, rcpt = c(1, 2))
ok("gspuff", mass = 5, u = 2, h = 10, rcpt = arc, t = 300, stab = "B", mixh = 200)
ac0 <- matrix(0, 7, 6)
ac0[3, 3] <- 1
ac0[4, 4] <- 2
asrc <- matrix(0, 7, 6)
asrc[5, 3] <- 0.5
ok("advdif2d", c0 = ac0, u = 0.7, v = -0.4, kx = 0.2, ky = 0.1, dx = 1, dy = 1, dt = 0.5, nsteps = 6, source = asrc)
ok("advdif2d", c0 = ac0, u = 2, v = 1, kx = 0.5, ky = 0.5, dx = 1, dy = 1, dt = 1, nsteps = 3)
bad("advdif2d", c0 = ac0, u = 0, v = 0, kx = 0, ky = 0, dx = 1, dy = 1, dt = 1, nsteps = 1e9)
ok("lagrprt", n = 5, x0 = 1, y0 = 2, u = 0.8, v = 0.1, kx = 0.5, ky = 0.3, dt = 0.7, nsteps = 4, seed = 5,
   grid = c(-2, 6, -2, 6, 4, 4))
ok("lagrprt", n = 3, x0 = 0, y0 = 0, u = 1, v = 1, kx = 1, ky = 1, dt = 1, nsteps = 2, seed = -7)
bad("lagrprt", n = 0, x0 = 0, y0 = 0, u = 1, v = 1, kx = 1, ky = 1, dt = 1, nsteps = 1)
ok("plmfld", q = 100, u = 5, h = 50, x = c(100, 600, 2000), y = c(-100, 0, 100))
bad("plmfld", q = 100, u = 5, h = 50, x = c(3, 1), y = 0)
ok("crfpm25", expo = c(5, 12, 20))
ok("crfpm25", expo = 17, y = "stroke", refconc = 5)
bad("crfpm25", expo = 10, y = "nope")
ok("crfno2", expo = c(15, 40), y = "chiast")
ok("crfno2", expo = 33, betaper10 = -0.02)
bad("crfno2", expo = 20, y = "nope")
ok("attrfrac", rr = 1.08, expoprev = 0.7, rrlo = 1.06, rrup = 1.09)
bad("attrfrac", rr = 1.5, expoprev = 1.2)
ok("mortdisp", expdel = 5, pop = 2.8e6, brate = 0.0075, betaperunit = log(1.08) / 10)
ok("mortdisp", expdel = -3, pop = 1e5, brate = 5e-6, betaperunit = 0.01, rare = TRUE)
bad("mortdisp", expdel = 1, pop = 1, brate = 0.5, betaperunit = 0.001)
ok("polbrdn", expomean = 12, expoprev = 0.9, brate = 0.005, pop = 5e5, poll = "PM2.5", y = "ihd")
ok("polbrdn", expomean = 30, expoprev = 1, brate = 0.012, pop = 5e5, y = "chiast", refconc = 8)
bad("polbrdn", expomean = 4, expoprev = 1, brate = 0.008, pop = 1e6, poll = "O3")
ok("polbrdarea", areatab = data.frame(fsa = c("A", "B", "C"), exposure = c(18, 28, 9),
                                      population = c(60000, 40000, 5000), brate = c(0.008, 0.007, 0.01)))
bad("polbrdarea", areatab = data.frame(fsa = "A", exposure = 20, population = 1, brate = 800))
ok("expconcidx", data = data.frame(h = c(30, 25, 25, 15, 10, 12), i = c(1, 2, 2, 3, 5, 4)), expo = "h",
   income = "i")
bad("expconcidx", data = data.frame(e = 1:2, i = 1), expo = "e", income = "i")

# One-pass summaries
skx <- c(3.5, -1, 2.25, 8, 1e8 + 0.5, 1e8 - 0.25, NA, 4, 0.125, 7)
ska <- smryupd(onlsmry(skx[1:4]), skx[5:7])
ok("onlsmry", x = skx[1:4])
ok("onlsmry")
bad("onlsmry", x = c("a", "b"))
ok("smryupd", acc = onlsmry(skx[1:4]), x = skx[5:7])
bad("smryupd", acc = list(), x = 1)
ok("smrymrg", a = ska, b = onlsmry(skx[8:10]))
ok("smrystat", acc = smrymrg(ska, onlsmry(skx[8:10])))
ok("smrystat", acc = onlsmry(c(1, 2)))
ok("rsvidx", n = 1000, k = 7, seed = 1)
ok("rsvidx", n = 50, k = 5, seed = 123456789)
bad("rsvidx", n = 10, k = -1)
ok("rsvsmp", x = letters, k = 4, seed = 2)
ok("dstnskch", x = c("a", "b", "c", "a", "1", 1, 2.5, NA), p = 4)
ok("dstnskch", x = c("d", "e"), p = 4, regstr = dstnskch(c("a", "b"), p = 4))
bad("dstnskch", x = 1, p = 3)
ok("dstncnt", regstr = dstnskch(as.character(1:3000), p = 6))
bad("dstncnt", regstr = c(1, 99))
ok("skchmrg", a = dstnskch(1:5, p = 4), b = dstnskch(6:9, p = 4))
bad("skchmrg", a = dstnskch(1:10, p = 4), b = dstnskch(1:10, p = 5))
# First Nations, Metis and Inuit health and wellbeing, published tables,
# the census profile and the Calls to Action (made-up numbers)
igd <- data.frame(g = c("First Nations", "Metis", "Total Indigenous", "Inuit", "Non-Indigenous"),
                  est = c(45.8, 52.8, 49, 47.8, 61.5), lo = c(43.9, 50.6, 47.5, 43.1, 61),
                  hi = c(47.7, 55, 50.5, 52.5, 62), q = c("", "", "", "E", ""))
ok("sctidy", x = data.frame(GEO = "Canada", VALUE = c("45.8", "12", "", "3"), STATUS = c("", "E", "x", "F")))
ok("sctidy", x = data.frame(GEO = "Canada", Coordinate = c("1.1", "1.2"), `Age (2):Total[1]` = c(10, 20),
                           Symbol = c("", "x"), check.names = FALSE))
bad("sctidy", x = data.frame(a = 1))
ok("qltyrcd", x = c("", "np", "*", NA), suppr = "np", caution = "*")
ok("qltyrcd", x = c("E", "Q"))
bad("qltyrcd", x = c(1, 2))
ok("limfrom", est = c(24.1, 31.5), cv = c(8, 12))
ok("limfrom", est = 1520, moe = 210, conflvl = 0.9)
bad("limfrom", est = 1)
ok("indgap", data = igd, grp = "g", est = "est", lower = "lo", upper = "hi", qual = "q")
ok("indgap", data = igd, grp = "g", est = "est", lower = "lo", upper = "hi", pooled = "Total Indigenous",
   ref = "Metis")
bad("indgap", data = igd, grp = "zz", est = "est", lower = "lo", upper = "hi")
ok("genexpo", parent = c(TRUE, FALSE, TRUE, NA), gparent = c(TRUE, FALSE, FALSE, TRUE))
ok("genexpo", parent = c(TRUE, NA), gparent = c(NA, TRUE), partial = "lowerbnd")
bad("genexpo", parent = 1, gparent = TRUE)
ok("exptrnd", cases = c(60, 95, 70), n = c(500, 520, 300))
ok("exptrnd", cases = c(0, 5), n = c(10, 10), dose = c(1, 3), conflvl = 0.9)
bad("exptrnd", cases = 1, n = 10)
ok("trccmpl")
ok("trcprog")
ok("trcprog", done = c(2, 6, 9, 13), yearend = 2016:2019, total = 20)
bad("trcprog", done = 1, yearend = 2010)
iwm <- data.frame(community = rep(c("Lakeside", "Riverbend"), each = 2), item = rep(c("Hopeful", "Meaning"), 2),
                  pct = c(64, 70, 58, 75), lo = c(59, 65, 52, 70), hi = c(69, 75, 64, 80))
ok("mhcont", data = iwm, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi",
   hope = "Hopeful", meaning = "Meaning", purp = "Meaning")
bad("mhcont", data = iwm, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi")
ok("cwbscr", incpercap = c(18000, 31000), highsch20p = c(58, 81), uni25p = c(6, 19), notcro = c(71, 95),
   nomajrep = c(62, 90), partic = c(55, 68), empl = c(78, 92))
bad("cwbscr", incpercap = -1, highsch20p = 58, uni25p = 6, notcro = 71, nomajrep = 62, partic = 55, empl = 78)
icn <- data.frame(region = c("North", "South", "East", "West", "Coast", "Interior"), language = c(41, 18, 27, 33, 22, 38),
                  wellbeing = c(68, 55, 60, 66, 57, 63))
icn$llo <- icn$language - 3
icn$lhi <- icn$language + 3
icn$wlo <- icn$wellbeing - 2
icn$whi <- icn$wellbeing + 2
ok("connassoc", data = icn, place = "region", x = "language", y = "wellbeing", reps = 200)
ok("connassoc", data = icn, place = "region", x = "language", y = "wellbeing", xlo = "llo", xup = "lhi",
   ylo = "wlo", yup = "whi", reps = 200, seed = 3)
bad("connassoc", data = icn, place = "region", x = "language", y = "wellbeing", xlo = "llo")
ihc <- data.frame(community = "Lakeside", item = rep(c("Good mental health", "Anxiety"), 2),
                  survey = rep(c("2018", "2023"), each = 2), pct = c(55, 14, 61, 11), lo = c(50, 11, 56, 8),
                  hi = c(60, 17, 66, 14))
ok("healchg", data = ihc, grp = "community", ind = "item", prd = "survey", est = "pct", lower = "lo", upper = "hi",
   from = "2018", to = "2023", well = "Good mental health", distress = "Anxiety")
bad("healchg", data = ihc, grp = "community", ind = "item", prd = "survey", est = "pct", lower = "lo",
    upper = "hi", from = "2019", to = "2023")
ici <- data.frame(community = rep(c("Lakeside", "Riverbend"), each = 3),
                  item = rep(c("Speaks the language", "Has a doctor", "Food insecure"), 2),
                  pct = c(38, 72, 21, 25, 64, 30), lo = c(33, 67, 17, 20, 58, 25), hi = c(43, 77, 25, 30, 70, 35))
ok("commidx", data = ici, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi",
   w = c("Speaks the language" = 2, "Has a doctor" = 1, "Food insecure" = -1))
bad("commidx", data = ici, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi", w = c(1, 2))
iwp <- data.frame(community = rep(c("Lakeside", "Riverbend"), each = 4),
                  item = rep(c("Hopeful", "Satisfied", "Rarely hopeful", "Low satisfaction"), 2),
                  pct = c(64, 71, 12, 18, 58, 66, 16, 22), lo = c(59, 66, 9, 14, 52, 60, 12, 17),
                  hi = c(69, 76, 15, 22, 64, 72, 20, 27))
ok("welprof", data = iwp, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi",
   well = c("Hopeful", "Satisfied"), distress = c("Rarely hopeful", "Low satisfaction"),
   domain = c(Hopeful = "spirit"))
bad("welprof", data = iwp, grp = "community", ind = "item", est = "pct", lower = "lo", upper = "hi",
    well = "Hopeful", distress = "Hopeful")
iw <- welprof(iwp, "community", "item", "pct", "lo", "hi", well = c("Hopeful", "Satisfied"),
              distress = c("Rarely hopeful", "Low satisfaction"))
ok("frmaud", x = iw)
bad("frmaud", x = data.frame(a = 1))
ok("dualcont", x = iw, well = "Hopeful", distress = "Rarely hopeful")
bad("dualcont", x = iw, well = "Rarely hopeful", distress = "Hopeful")
ok("scests", x = data.frame(GEO = "Canada", Characteristics = c("Rate", "Coefficient of variation for rate"),
                           VALUE = c(12.5, 8), STATUS = ""))
ok("scests", x = data.frame(GEO = "Canada", g = rep(c("A", "B"), each = 3),
                           Statistics = rep(c("Percent", "Low 95% confidence interval",
                                              "High 95% confidence interval"), 2),
                           VALUE = c(47.1, 45.2, 49, 50.3, 48.1, 52.5), STATUS = ""))
bad("scests", x = data.frame(GEO = "Canada", Statistics = c("Percent", "Percent"), VALUE = c(1, 2)))
ioc <- data.frame(group = rep(c("Nation A", "Everyone else"), each = 4),
                  rating = rep(c("Poor", "Fair", "Good", "Excellent"), 2), pct = c(10, 25, 40, 25, 6, 20, 44, 30),
                  lo = c(7, 21, 35, 21, 5, 18, 42, 28), hi = c(13, 29, 45, 29, 7, 22, 46, 32))
ok("ordcmp", data = ioc, grp = "group", cat = "rating", est = "pct", lower = "lo", upper = "hi",
   lvls = c("Poor", "Fair", "Good", "Excellent"), ref = "Everyone else")
bad("ordcmp", data = ioc, grp = "group", cat = "rating", est = "pct", lower = "lo", upper = "hi",
    lvls = "Poor", ref = "Everyone else")
ok("stdprev", est = c(12, 9, 6, 8, 6, 4), lower = c(9, 6, 3, 7, 5, 3), upper = c(15, 12, 9, 9, 7, 5),
   std = c(500, 300, 200), grp = rep(c("Nation A", "Everyone else"), each = 3))
bad("stdprev", est = 1, lower = 0, upper = 2, std = 0)
ihb <- data.frame(community = rep(c("Lakeside", "Riverbend"), each = 4),
                  domain = rep(c("body", "mind", "spirit", "family"), 2), pct = c(62, 58, 70, 75, 55, 61, 66, 80),
                  lo = c(57, 53, 65, 70, 49, 55, 60, 75), hi = c(67, 63, 75, 80, 61, 67, 72, 85))
ok("holbal", data = ihb, grp = "community", domain = "domain", est = "pct", lower = "lo", upper = "hi", reps = 200)
ok("holbal", data = ihb, grp = "community", domain = "domain", est = "pct", lower = "lo", upper = "hi", reps = 100,
   seed = 4, w = c(body = 1, mind = 1, spirit = 2, family = 1), conflvl = 0.9)
bad("holbal", data = ihb, grp = "community", domain = "domain", est = "pct", lower = "lo", upper = "hi", reps = 1)
ipp <- data.frame(GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
                  AGE_GROUP = c(rep("Total - Age groups", 3), rep("25 years and over", 3), rep("25 to 64 years", 2)),
                  RESIDENCE = "Total - Residence on or off reserve", CHARACTERISTIC_ID = 1:8,
                  CHARACTERISTIC_NAME = c(
                    "Total - Age groups of the population - 100% data",
                    paste("Number of total income recipients aged 15 years and over in private households in",
                          "2020 - 25% sample data"),
                    "Average total income in 2020 among recipients ($)",
                    paste("Total - Highest certificate, diploma or degree for the population aged 25 years and",
                          "over in private households - 25% sample data"),
                    "  No certificate, diploma or degree", "  Bachelor's degree or higher", "Participation rate",
                    "Unemployment rate"),
                  `IND_IDENT (2)` = c(1200, 900, 38000, 700, 150, 120, 62, 14), SYMBOL_IND_IDENT = c("", "", "", "",
                    "E", "", "", ""),
                  NON_IND = c(5000, 4000, 52000, 3500, 400, 900, 66, 6), SYMBOL_NON_IND = "", check.names = FALSE)
ok("cnsprof", x = ipp)
bad("cnsprof", x = data.frame(a = 1))
ip <- cnsprof(ipp)
ok("profshr", x = ip, num = "No certificate, diploma or degree")
ok("profshr", x = ip, num = 6, den = 4, smpfrac = 0.2, roubase = 10)
bad("profshr", x = ip, num = 99)
ist <- data.frame(GEO_NAME = "Lakeside", DGUID = "X0001", GENDER = "Total - Gender",
                  AGE_GROUP = c("Total - Age groups", "0 to 14 years", "15 to 64 years", "65 years and over"),
                  RESIDENCE = "Total - Residence on or off reserve", CHARACTERISTIC_ID = 1,
                  CHARACTERISTIC_NAME = "Total - Age groups of the population - 100% data",
                  `IND_IDENT (2)` = c(1200, 400, 700, 100), SYMBOL_IND_IDENT = "", NON_IND = c(5000, 800, 3300, 900),
                  SYMBOL_NON_IND = "", check.names = FALSE)
ok("cnsstd", x = cnsprof(ist), geo = "Lakeside", ages = c("0 to 14 years", "15 to 64 years", "65 years and over"))
bad("cnsstd", x = cnsprof(ist), geo = "Riverbend")
ok("cnscwb", x = ip, groups = c("Indigenous identity", "Non-Indigenous"), minpop = 0)
bad("cnscwb", x = data.frame(a = 1))
ok("trcind")
ilt <- data.frame(period = "2017/2019", group = c("First Nations", "Inuit", "Non-Indigenous"),
                  estimate = c(74.2, 70.9, 81.8), lower = c(73.6, 69.1, 81.6), upper = c(74.8, 72.7, 82.0),
                  quality = "ok")
ok("trchlthgap", lifeexp = ilt)
# a numeric period: copied to every row of the result
ok("trchlthgap", lifeexp = transform(ilt, period = 2019))
bad("trchlthgap")
ok("polex", seed = 3)
bad("polex", seed = 1.5)

# Geostatistics
gp <- cbind(x = c(0, 1, 0, 1, 2, 0.5, 1.5, 2.2), y = c(0, 0, 1, 1, 2, 1.7, 0.4, 1.1))
gv <- c(1, 2, 3, 5, 4, 2.5, 3.1, 4.4)
ok("vgemp", x = gp, v = gv)
ok("vgemp", x = gp, v = gv, nbins = 4, cutoff = 2)
bad("vgemp", x = gp, v = gv, nbins = 1e9)
bad("vgemp", x = gp, v = gv[1:3])
ok("vgfit", x = gp, v = gv)
ok("vgfit", x = gp, v = gv, model = "spherical")
bad("vgfit", x = gp, v = gv, model = "linear")
ok("okrge", x = gp, v = gv, at = rbind(c(0.5, 0.5), c(1.8, 1.9)))
ok("okrge", x = gp, v = gv, at = rbind(c(0.5, 0.5)),
   vgm = list(model = "gaussian", nugget = 0.1, psill = 1, range = 1.2))
bad("okrge", x = gp, v = gv, at = rbind(c(0.5, 0.5)),
    vgm = list(model = "exponential", nugget = 0, psill = 1, range = 0))
ok("bkrig", x = gp, v = gv, blocks = list(rbind(c(0, 0), c(1, 1)), rbind(c(1.2, 0.2), c(1.4, 0.3), c(1.6, 0.1))),
   nq = 9, range = 0.8)
bad("bkrig", x = gp, v = gv, blocks = list(1:3))
bad("bkrig", x = gp, v = gv, blocks = list(rbind(c(0, 0), c(1, 1))), nq = 2147483647)
ok("krigcv", x = gp, v = gv, nugget = 0.1, sill = 1.2, range = 0.7)
bad("krigcv", x = gp, v = gv, nugget = 2, sill = 1)
ok("trsurf", x = gp, v = gv)
ok("trsurf", x = gp, v = gv, order = 1)
bad("trsurf", x = gp, v = gv, order = 5)
ok("sglmx", v = gv, z = cbind(1, gp[, 1]), x = gp)
ok("sglmx", v = gv, x = gp)
bad("sglmx", v = gv, z = cbind(1, 1:3), x = gp)
ok("smixm", v = gv, z = cbind(1, gp[, 1]), x = gp)
bad("smixm", v = gv, z = "a", x = gp)
gs <- spobj(data.frame(lon = c(-79.4, -79.3, -79.5), lat = c(43.6, 43.7, 43.65), v = c(1, 2, 3)))
ok("spobj", data = data.frame(lon = c(-79.4, -79.3), lat = c(43.6, 43.7), v = c(1, 2)))
ok("spobj", data = data.frame(e = c(1, 2), n = c(3, 4)), coords = c("e", "n"), crs = 3857)
bad("spobj", data = data.frame(lon = 1, lat = 1), crs = 1.5)
ok("spxy", x = gs)
bad("spxy", x = list(a = 1))
ok("crsof", x = gs)
ok("crsof", x = spobj(data.frame(e = 1, n = 2), c("e", "n"), crs = 2958))
ok("lfspec", x = gs)
ok("lfspec", x = spobj(data.frame(e = c(-8800000, -8810000), n = c(5400000, 5410000)), c("e", "n"), crs = 3857))
bad("lfspec", x = spobj(data.frame(e = 1, n = 2), c("e", "n"), crs = 2958))
ok("spsamp", x = gs, size = 2, seed = 1)
ok("spsamp", x = gs, size = 3, dens = FALSE, seed = 2)
bad("spsamp", x = gs, size = 4, seed = 1)

# Spatial epidemiology
ok("bayout", x = c(3, 5, 2, 4, 6, 3, 4, 12), w = 6, at = 8)
ok("bayout", x = c(4, 5, 3, NA, 6, 2, 7, 5, 9, 4, 3, 15), freq = 4, b = 1, w = 1, acty = FALSE, alpha = 0.1)
bad("bayout", x = c(1, 2, 3), w = 6)
er <- cbind(x = c(0, 1, 5), y = c(0, 1, 5))
es <- cbind(x = c(0.5, 4), y = c(0.5, 4))
ok("bufexp", x = er, src = es, r = 1.5)
ok("bufexp", x = cbind(lon = c(-79.4, -79.38), lat = c(43.65, 43.7)),
   src = cbind(lon = c(-79.39, -79.5), lat = c(43.66, 43.6)), r = 5000, w = c(2, 1), lonlat = TRUE)
bad("bufexp", x = er, src = es, r = -1)
ok("kerexp", x = er, src = es, h = 2)
ok("kerexp", x = er, src = es, h = 2, w = c(1, 3), kern = "quartic")
bad("kerexp", x = er, src = es, h = 2, kern = "epan")
ew <- rbind(c(0, 1, 1, 0), c(1, 0, 1, 0), c(1, 1, 0, 1), c(0, 0, 1, 0))
ok("lerprc", w = ew, rho = 0.5, tau = 2)
bad("lerprc", w = ew, rho = 2)
ok("bym2sc", w = ew, tau = 2, phi = 0.4)
bad("bym2sc", w = rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 0)))
ey <- c(2, 0, 5, 3, 0, 7, 1, 4)
ee <- c(2.5, 1, 4, 3, 1.5, 5, 1.2, 3.3)
ex <- c(0.1, -0.3, 0.5, 0.2, -0.6, 0.9, -0.1, 0.4)
ok("ecolrg", y = ey, e = ee, x = ex)
ok("ecolrg", y = ey, e = ee, x = ex, family = "negbin")
ok("ecolrg", y = ey, e = ee, x = ex, family = "zip", maxit = 50, tol = 1e-10)
bad("ecolrg", y = ey, e = ee, family = "gamma")
ok("spbias", d = rbind(c(0, 1, 2), c(1, 0, -1), c(2, -1, 0)), b = c(1, 0.5, 0.25))
bad("spbias", d = rbind(c(0, 1.5)), b = 1)
ok("spsir", i0 = c(10, 0, 0), w = rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0)), pop = c(1000, 800, 1200),
   beta = 0.6, tmax = 4, dt = 1)
ok("spsir", i0 = c(5, 1), tmax = 2, dt = 0.5)
bad("spsir", i0 = c(5, 1), pop = c(3, 3))
fi <- 0:11
ft <- 0.5 + ((fi * 7) %% 11) / 5
fw <- rbind(c(0, 1, 0), c(1, 0, 1), c(0, 1, 0))
ok("sfrail", time = ft, event = as.numeric(fi %% 4 != 0), x = sin(fi), reg = fi %% 3 + 1, w = fw)
ok("sfcure", time = ifelse(fi %% 4 == 0, 5, ft), event = as.numeric(fi %% 4 != 0 & fi %% 5 != 1), reg = fi %% 3 + 1,
   w = fw, tau = 2)
bad("sfrail", time = ft, event = fi, reg = fi %% 3 + 1, w = fw)
ok("evfilt", v = sin(1:8 / 2), w = 1 * (abs(outer(1:8, 1:8, "-")) == 1))
ok("evfilt", v = sin(1:8 / 2), w = 1 * (abs(outer(1:8, 1:8, "-")) == 1), x = cos(1:8), method = "aic", cand = "all")
bad("evfilt", v = 1:2, w = diag(2))
# Summary statistics in compiled form
px <- c(2, 4, 4, 4, 5, 5, 7, 9, 1.5, 3.25)
py <- c(1, 3, 2, 5, 4, 6, 8, 7, 9, 10)
ok("cmn", x = px)
bad("cmn", x = c("a", "b"))
ok("cvari", x = px)
ok("ccorr", x = px, y = py)
bad("ccorr", x = 1:3, y = 1:4)
ok("csdev", x = px, ddof = 0)
bad("csdev", x = px, ddof = -1)
ok("cdistm", a = px, b = py)
ok("cnormpdf", x = c(-1, 0, 2.5), mean = 0.3, sd = 1.7)
bad("cnormpdf", x = 0, sd = 0)
ok("cnormlpdf", x = c(-1, 0, 50), mean = 0.3, sd = 1.7)
ok("cmom", x = px)
ok("cmom", x = c(1, 2))
ok("cqnt", x = px, probs = c(0, 0.1, 0.33, 0.5, 0.9, 1))
ok("cqnt", x = px)
bad("cqnt", x = px, probs = 1.5)
ok("cmed", x = px)
ok("cmad", x = px)
ok("cmad", x = px, const = 1)
bad("cmad", x = px, const = -1)
ok("ciqr", x = px)
ok("ctkyfnc", x = px, k = 2)
ok("ctrmmean", x = c(px, 100), trim = 0.2)
bad("ctrmmean", x = px, trim = 0.7)
ok("cwinsmean", x = c(px, 100), trim = 0.2)
ok("cwtd", x = px, w = py)
bad("cwtd", x = px, w = -py)
ok("cspr", x = px, y = py)
ok("cmidrk", x = c(3, 1, 3, 2))
ok("ccovm", x = cbind(px, py))
ok("ccovm", x = data.frame(a = px, b = py))
bad("ccovm", x = px)
ok("cbootmean", x = px, B = 6, seed = 7)
bad("cbootmean", x = px, B = 0)
ok("cipw", treat = c(1, 0, 1, 0), ps = c(0.5, 0.25, 0.002, 0.995))
ok("cipw", treat = c(TRUE, FALSE), ps = c(0.4, 0.4), trimlo = 0.1, trimhi = 0.9)
bad("cipw", treat = c(1, 2), ps = c(0.5, 0.5))
ok("cgamcdf", shape = 3.5, x = c(0.5, 1, 4, 12, 0, -1))
bad("cgamcdf", shape = 0, x = 1)
ok("cunif", n = 4, seed = 42)
bad("cunif", n = 2, seed = 2^53)

# Text to plain ASCII
ok("toascii", x = c("Prof. \u00c1ngela Zorro Medina", "Stra\u00dfe \u00c6sir \u0152uvre \u00de\u00fe",
                     "\u0141\u00f3d\u017a \u0160koda", NA, "a\tb", "\u00d7\u00ff"))
ok("toascii", x = c(1.5, 2))
ok("asciifb", x = c("\u00c1ngela", NA, "plain"))
ok("asciifb", x = c("\u00c1ngela", NA), force = TRUE)
bad("wrtxtfb", text = "x", path = c("a", "b"))

# Shares with Wilson intervals
ss <- data.frame(division = c("North", "South", "East", "North", "East"), yr = c(1, 1, 1, 2, 2),
                 stops = c(412, 77, 3, 20, NA))
ok("shrs", x = ss, count = "stops", by = c("division", "yr"))
ok("shrs", x = ss, count = "stops", total = 10000)
ok("shrs", x = ss, count = "stops", total = "stops", conflvl = 0.8)
ok("shrs", x = c(a = 0, b = 7, c = 3), conflvl = 0.9)
ok("shrs", x = c(3, 7))
bad("shrs", x = c(5, 2), total = 4)
bad("shrs", x = ss, count = "stops", by = "zz")

# Categorical integrity
cm <- c("1" = "White", "2" = "Black", "3" = "Indigenous")
ok("deccodes", x = c(3, 1, NA, 2, 2), vlabs = cm)
bad("deccodes", x = c(1, 4), vlabs = cm)
bad("deccodes", x = c(1, NA), vlabs = cm, keepna = FALSE)
ok("grdrcd", x = c("M", "F", "F", "X", NA), map = c(M = "Male", F = "Female"), keep = "X")
bad("grdrcd", x = c("M", "U"), map = c(M = "Male"))
ok("relbl", x = c("b", "a", "c", "a"), map = c(c = "C", a = "A", b = "B"))
ok("relbl", x = c(2, 10, 1), map = c("1" = "one", "2" = "two", "10" = "ten"))
bad("relbl", x = c("W", "Q"), map = c(W = "White"))
ok("grdlvls", x = c("b", "a", NA), lvls = c("b", "a"), ref = "b")
bad("grdlvls", x = "c", lvls = c("b", "a"))
bad("grdlvls", x = "a", lvls = c("b", "a"), ref = "a")
ok("declbld", x = c(10, 2, 2, 10), vlabs = c("10" = "Ten", "2" = "Two"))
bad("declbld", x = c(1, 2))
ok("grdbin", x = c(0, 1, NA), col = "d")
bad("grdbin", x = c("treated", "control"), col = "d")
bad("grdbin", x = c(0, 2, 0.5), col = "d")
ok("vrfrcd", orig = c("M", "F", "F", NA, "X"), recd = c("Male", "Female", "Female", NA, "X"),
   decl = c(M = "Male", F = "Female"))
bad("vrfrcd", orig = c("M", "F"), recd = c("Female", "Male"), decl = c(M = "Male", F = "Female"))
bad("vrfrcd", orig = c("M", "M"), recd = c("Male", "Female"), decl = c(M = "Male"))
cx <- rep(c("White", "Black", "Indigenous"), c(5, 3, 2))
ok("vrfmrg", x = cx, pub = c(White = 3, Black = 5, Indigenous = 2), strict = FALSE)
ok("vrfmrg", x = cx, pub = c(White = 5, Black = 3, Indigenous = 2))
ok("vrfmrg", x = c(cx, "Other"), pub = c(White = 5, Black = 3, Indigenous = 2), strict = FALSE)
bad("vrfmrg", x = cx, pub = c(White = 3, Black = 5, Indigenous = 2))
ok("xfrvrf", imported = c(1, 2, 2, 1, 1), soucnts = c(No = 3, Yes = 2), codebook = c("1" = "No", "2" = "Yes"))
ok("xfrvrf", imported = c("No", "Yes", "Yes"), soucnts = c(No = 2, Yes = 1), strict = FALSE)
bad("xfrvrf", imported = c("No", "Yes", "Yes"), soucnts = c(No = 2, Yes = 1))

# Region-map checks
rcw <- data.frame(inst = c("North Jail", "South Jail", "East Jail", "North Jail", "West Jail"),
                  cd = c("3557", "3520", "2406", "3560", ""))
ok("rgnmapintg", map = rcw, unit = "inst", reg = "cd")
ok("rgnmapintg", map = rcw, unit = "inst", reg = "cd", regs = c("3557", "3520", "3506"))
bad("rgnmapintg", map = rcw, unit = "facility", reg = "cd")
rpub <- data.frame(inst = c("North Jail", "South Jail", "Hill Jail"), cd = c("3557", "3520", "3506"),
                   lat = c(44.39, 43.65, 45))
robs <- data.frame(inst = c("North Jail", "South Jail", "Far Jail"), cd = c("3557", "3521", "3506"),
                   lat = c(44.39 + 1e-12, 43.7, 45))
ok("rgnmapcmp", pub = rpub, obs = robs, unit = "inst")
ok("rgnmapcmp", pub = rpub, obs = robs, unit = "inst", cols = "cd")
bad("rgnmapcmp", pub = rpub, obs = robs, unit = "inst", cols = "x")
ok("rgnmapalt", map = rpub, unit = "inst", reg = "cd",
   route = c("North Jail" = "3557", "South Jail" = "3521", "Hill Jail" = "3519"), known = "Hill Jail")
ok("rgnmapalt", map = rpub, unit = "inst", reg = "cd",
   route = data.frame(inst = c("Hill Jail", "South Jail"), cd = c("3506", "9")))
bad("rgnmapalt", map = rpub, unit = "inst", reg = "cd", route = c("3557", "3520"))

# The drift battery
dx <- c(3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5)
dy <- c(2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5)
ok("drftks", x = dx, y = dy)
ok("drftks", x = 1:300, y = (1:300) + 40)
bad("drftks", x = c(1, Inf), y = 1)
ok("drftpsi", x = dx, y = dy)
ok("drftpsi", x = dx, y = dy, bins = 4, eps = 0.01)
ok("drftpsi", x = c(1, 1, 1), y = c(9, 9, 9))
bad("drftpsi", x = dx, y = dy, bins = 1)
ok("drftchisq", obs = c(a = 40, b = 60, c = 100), expd = c(a = 50, b = 30, c = 20))
ok("drftchisq", obs = c("a", "a", "b", "b"), expd = c("a", "b"))
ok("drftchisq", obs = c(a = 5), expd = c(a = 2))
bad("drftchisq", obs = c(1, 2), expd = c(a = 1))
ok("drfthomog", x = c(a = 30, b = 20), y = c(a = 25, b = 25))
ok("drfthomog", x = c(a = 1, b = 1), y = c(c = 1, d = 1))
bad("drfthomog", x = c(a = 0), y = c(a = 3))

# Change tables out
yt <- yoypx(c(812, 845, 790, 120, 131, 98), rep(2021:2023, 2), rep(c("Assault", "Robbery"), each = 3))
ok("fylbl", endyear = c(2023, 2001, NA))
ok("fylbl", endyear = "2023", sep = "-", short = FALSE)
bad("fylbl", endyear = 2023.5)
ok("fylbl", endyear = c(1e300, Inf, NaN, 2023))
ok("fylbl", endyear = "3000000000")
ok("yoylbl", x = yt, labels = c("2021" = "first"))
ok("yoylbl", x = yt, labels = c("a", "b", "c"))
bad("yoylbl", x = yt, labels = c("a", "b"))
ok("yoysmry", x = yt)
ok("yoysmry", x = yt, lag = 2)
bad("yoysmry", x = data.frame(a = 1))
ok("yoycsv", x = yt)
ok("yoycsv", x = data.frame(code = c("046", "", "a,b"), n = c(1 / 3, NA, 1e-20)), meta = "simulated", na = "NA")
ok("yoytsv", x = yt, digits = 2)
bad("yoytsv", x = yt, file = strrep("x", 5000))
ok("yoyjs", x = yt)
ok("yoyjs", x = data.frame(code = c("046", "q\"x"), ok = c(TRUE, NA)), pretty = FALSE, meta = "m")
bad("yoywr", x = data.frame(n = 1), file = "x.out", format = "nope")
ok("yoymd", x = yt)
ok("yoymd", x = yoypx(c(5, 6.25), 1:2, units = "percent"), digits = 2, align = FALSE, meta = c("a", "b"))
bad("yoymd", x = yt, digits = 20)

# The footprint of a computation (co2fp() evaluates an R expression: R only)
ok("co2int", loc = "CA-ON")
ok("co2int", loc = "CA-XX")
ok("co2int", loc = 120)
bad("co2int", loc = "Atlantis")
ok("co2inttab")
ok("dtctloc", envloc = "", codloc = "", tz = "Europe/Monaco")
ok("dtctloc", envloc = "", codloc = "", tz = "", locale = "en_IN.UTF-8")
ok("dtctloc", envloc = "CA-QC", tz = "Europe/Paris")
ok("dtctloc", envloc = "Atlantis", codloc = "", tz = "Antarctica/Troll")
bad("dtctloc", tz = c("a", "b"))
ok("raplok", rapldir = "/nonexistent/powercap")
ok("fpeqv", co2eg = 1200)
bad("fpeqv", co2eg = -1)

# The exposure-response estimate (simulated areas)
set.seed(1)
exd <- data.frame(x1 = round(rnorm(40), 3), x2 = round(rnorm(40), 3))
exd$exposure <- round(20 + 2 * exd$x1 + rnorm(40), 3)
exd$asthma <- round(5 + 0.3 * exd$exposure + exd$x2 + rnorm(40), 3)
exd2 <- exd
exd2$x2[c(3, 7)] <- NA
exd2$asthma[11] <- NA
ok("exprspplr", data = exd[1:30, ], y = "asthma", expo = "exposure", confs = c("x1", "x2"), nboot = 12)
ok("exprspplr", data = exd2, y = "asthma", expo = "exposure", confs = "x1", nboot = 10, rstate = 2^35 + 1)
bad("exprspplr", data = exd[1:12, ], y = "asthma", expo = "exposure", confs = "x1", nboot = 5)
bad("exprspplr", data = exd[1:8, ], y = "asthma", expo = "exposure", confs = "x1")
bad("exprspplr", data = exd[1:12, ], y = "asthma", expo = "exposure", confs = "zz")

# Data description
ok("inlhist", x = c(1, 2, 2, 3, 3, 3, 7.5, NA), bins = 5)
ok("freqtbl", data = data.frame(grade = c("a", "b", "b", "c", NA, "b")), col = "grade")
ok("freqtbl", data = c("c", "a", "b", "a", "d", "d", "d"), maxlvl = 2)
dsdc <- data.frame(a = c(1.5, 2, NA, 4, 5, 3.25), b = c(2, 1, 4, 3, 6, 5), cc = c(9, 7, 8, 6, 4, 5), g = c("x", "y", "x", "y", "x", "y"))
ok("cortbl", data = dsdc)
ok("cortbl", data = dsdc, meth = "pearson", minprs = 6)
ok("topcors", data = dsdc, n = 2)
dsde <- data.frame(keep = c(1, 2, NA), allna = c(NA, NA, NA), constant = c(7, 7, 7), s = c("x", NA, "y"))
ok("drpempty", data = dsde, quiet = TRUE)
ok("drpempty", data = data.frame(a = c(1, NA), b = c(2, NA)), which = "rows", quiet = TRUE)
ok("drpconst", data = dsde, quiet = TRUE)
ok("clncolnms", data = c("Total  Population (2021)", "% change", "dup", "dup", "R&D #", "camelCaseName", "2021 count"), quiet = TRUE)
ok("clncolnms", data = data.frame(`Total Pop` = 1:2, `per capita` = 3:4, check.names = FALSE), case = "upper_camel", quiet = TRUE)
dsmh <- data.frame(h = c(170, 165, 180, 155, 172, 168, 175, 160, 178, 169), w = c(70, 60, 80, 100, 72, 66, 75, 58, NA, 68))
ok("mahaoutl", data = dsmh)
ok("mahaoutl", data = dsmh, robust = FALSE, alpha = 0.05)
dsmr <- data.frame(outage = c(1, 2, NA, NA, NA, NA, 7, 8), scattered = c(1, NA, 3, 4, NA, 6, NA, NA))
ok("missruns", data = dsmr)
ok("missgrid", data = dsmr, height = 4, width = 3)
ok("missgrid", data = data.frame(a = c(1, NA), bb = c(NA, 2)))
dsmc <- data.frame(x = c(0.5, -1.2, 0.3, 1.8, -0.4, 0.9, -2.1, 0.0, 1.1, -0.7, 0.2, 1.4),
                 y = c(1.1, -0.8, NA, 2.2, 0.1, NA, -1.9, 0.4, NA, -0.2, 0.6, 1.9),
                 z = c(0.3, 0.7, -0.5, NA, 1.2, -1.1, 0.4, -0.3, 0.8, NA, -0.6, 0.1))
ok("mcartst", data = dsmc)
ok("misspat", data = data.frame(id = 1:6, a = c(NA, NA, 3:6), b = c(1:4, NA, NA)))
ok("duprows", data = data.frame(id = c(1, 2, 2, 3, 3, 3), value = c("a", "b", "b", "c", "d", "c")), quiet = TRUE)
ok("duprows", data = data.frame(id = c(1, 2, 2, 3, 3, 3), value = c("a", "b", "b", "c", "d", "c")), cols = "id", quiet = TRUE)
ok("misssmry", data = data.frame(a = c(1, NA, 3), b = c(NA, NA, 3), c = 1:3))
ok("profcols", data = data.frame(v = c(1.5, 2, 2, 8, NA, -3, 0), g = c("a", "b", "b", NA, "a", "b", "c")), qnts = c(0.1, 0.5))
dsown <- data.frame(group = c("Aboriginal", "Torres Strait Islander", "Both"), est = c(24.0, 31.0, 24.6), lo = c(22.1, 26.2, 22.8), hi = c(25.9, 35.8, 26.4))
ok("dstntst", data = dsown, grp = "group", est = "est", lower = "lo", upper = "hi", pooled = "Both")
bad("inlhist", x = c("a", "b"))
bad("inlhist", x = c(1, 2, 3), bins = 1e9)
bad("freqtbl", data = data.frame(g = c("a", "b")))
bad("cortbl", data = data.frame(a = 1:5, b = letters[1:5]))
bad("topcors", data = dsdc, minabs = 2)
bad("drpempty", data = c(1, NA))
bad("drpconst", data = dsde, quiet = "no")
bad("clncolnms", data = "a b", case = "kebab")
bad("mahaoutl", data = dsmh, alpha = 1)
bad("missruns", data = dsmr, minrun = 0)
bad("missgrid", data = dsmr, height = 0)
bad("mcartst", data = data.frame(x = 1:3, g = c("a", "b", "c")))
bad("misspat", data = data.frame())
bad("duprows", data = data.frame(id = 1:3), cols = "code", quiet = TRUE)
bad("misssmry", data = c(1, NA, 3))
bad("profcols", data = data.frame(v = 1:3), qnts = 2)
bad("dstntst", data = dsown, grp = "grp", est = "est", lower = "lo", upper = "hi")

# Distributions fitted by maximum likelihood (simulated)
set.seed(1)
dfg <- round(rgamma(60, shape = 2, rate = 0.5), 3)
set.seed(2)
dfc <- ifelse(runif(50) < 0.3, 0, rpois(50, 4))
dfc[c(3, 9)] <- NA
set.seed(3)
dfp <- round(rbeta(40, 2, 5), 4)
ok("dstfit", x = dfg[1:30], nboot = 5, seed = 3, fams = c("gamma", "weibull", "lognormal", "exponential"))
ok("dstfit", x = dfc, nboot = 4, fams = c("poisson", "negbin", "zip"))
ok("dstfit", x = dfp, nboot = 0, fams = c("beta", "normal", "pareto"))
bad("dstfit", x = c(-1, 2, 3, 4, 5), fams = "gamma")
bad("dstfit", x = c(1, 2, 3))

# Scaling and ideal points
set.seed(21)
sv_v <- matrix(rbinom(48, 1, plogis(outer(rnorm(8), rnorm(6), "*") * 3)), 8)
sv_v[2, 3] <- NA
sv_d <- as.matrix(dist(cbind(c(0, 1, 3, 4, 2), c(0, 2, 1, 3, 4))))
sv_z <- rbind(c(1, 3, 5, 6), c(2, 3, 4, 7), c(1, 2, 6, 6), c(3, 4, 5, 7), c(2, 4, 4, 6), c(1, 3, 6, 7))
ok("votagr", x = sv_v)
bad("votagr", x = "a")
ok("amscal", Z = sv_z)
bad("amscal", Z = matrix(1:3, 3))
ok("bbscal", X = cbind(sv_z, sv_z[, 2:1] + c(0.5, -0.5)), ndim = 1, minscale = 4)
bad("bbscal", X = matrix(1:4, 4))
ok("optcls", votes = sv_v, maxit = 5, nstart = 2)
bad("optcls", votes = sv_v, ndim = 0)
bad("optcls", votes = sv_v, ndim = 1e9)
ok("dblctr", D = sv_d)
bad("dblctr", D = matrix(1:6, 2))
ok("cmdsf", D = sv_d)
ok("smacf", D = sv_d, maxit = 20)
bad("smacf", D = sv_d, w = matrix(1, 2, 2))
bad("smacf", D = sv_d, ndim = 1e9)
ok("nmdsf", D = sv_d, maxit = 10)
bad("nmdsf", D = sv_d, ndim = 1e9)
ok("mdsfit", ev = c(4, 2, -1, 0.5))
ok("unfstr", xr = matrix(c(0, 1, 0, 1), 2), xs = matrix(c(1, 0, 2, 1, 1, 0), 3), D = matrix(c(1, 2, NA, 1, 1, 1), 2))
bad("unfstr", xr = matrix(0, 2, 2), xs = matrix(0, 3, 2), D = matrix(1, 3, 2))
ok("mlsmu", D = sv_z, maxit = 10, nstart = 2)
ok("smunf", D = sv_z, maxit = 10)
bad("mlsmu", D = sv_z, ndim = 1e9)
bad("smunf", D = sv_z, ndim = 1e9)
ok("ipsrec", xr = matrix(1:4, 2))
ok("nomutl", x = matrix(c(0, 0.5, -0.5, 0.2), 2), zy = matrix(c(0.3, 0.1), 1), zn = matrix(c(-0.3, 0), 1))
ok("nompr", x = 0.1, zy = 0.3, zn = -0.3)
bad("nompr", x = 0.1, zy = 0.3, zn = -0.3, beta = "a")
bad("nompr", x = c(0.1, 0), zy = 0.3, zn = -0.3)
ok("nomll", votes = matrix(c(1, 0, NA, 1), 2), x = matrix(c(0.2, -0.4)), zy = matrix(c(0.3, -0.1)), zn = matrix(c(-0.3, 0.4)))
ok("procrs", X = sv_d[, 1:2], target = sv_d[, 2:1])
bad("procrs", X = matrix(1:4, 2), target = matrix(1:6, 3))
ok("irtlik", votes = sv_v[1:3, 1:4], x = matrix(c(-1, 0, 1)), alpha = c(0, 0.5, -0.5, 0), beta = matrix(c(1, 2, -1, 0.5)))
ok("irtpst", chain = array(seq(0.1, 2.4, by = 0.1)^2, c(4, 3, 2)))
bad("irtpst", chain = 1:3)
ok("ordocl", Y = sv_z, ndim = 1, maxit = 5)
bad("ordocl", Y = sv_z, ndim = 1e9)
ok("anchv", y = c(1, 3, 5, NA, 2, 4), V = sv_z[, 1:3])
bad("anchv", y = 1:3, V = sv_z)
bad("anchv", y = c(1, 3, 5, NA, 2, 4), V = sv_z[, 1:3], ncat = 1e9)
ok("indscl", D = list(sv_d, sv_d * c(1, 1.2, 0.9, 1.1, 1)), maxit = 10)
bad("indscl", D = list(sv_d, sv_d[1:3, 1:3]))
bad("indscl", D = list(sv_d, sv_d), ndim = 1e9)
ok("nrmvec", x = sv_d[, 1:2], y = c(1, 2, NA, 4, 3))
ok("cutlns", normals = rbind(c(0.6, 0.8), c(1, 0), c(0, 0)), cuts = c(0.1, -0.2, 0))
bad("cutlns", normals = matrix(1, 2, 3), cuts = 1:2)
ok("dwnom", votes = sv_v, maxit = 5)
bad("dwnom", votes = sv_v, ndim = 1e9)
ok("nombt", votes = sv_v[1:4, 1:3], x = matrix(c(-1, -0.5, 0.5, 1)), normals = matrix(c(1, 1, -1)), cuts = c(0, 0.2, -0.1), nboot = 3)
ok("irtem", votes = sv_v, maxit = 20)
ok("irtem", votes = sv_v, ndim = 2, maxit = 5, conv = "abs")
bad("irtem", votes = sv_v, xprior = 1)
bad("irtem", votes = sv_v, ndim = 1e9)
ok("scboot", Z = sv_z, nboot = 3)
bad("scboot", Z = sv_z, method = "x")
bad("scboot", Z = sv_z, nboot = 1e9)
ok("wrdfsh", dtm = matrix(c(5, 3, 0, 2, 4, 1, 6, 2, 3, 1, 0, 7), 4), maxit = 10)
bad("wrdfsh", dtm = matrix(-1, 2, 2))
ok("bamscl", Z = sv_z, nsamp = 5, burn = 3)
bad("bamscl", Z = sv_z, polarity = 9)
bad("bamscl", Z = sv_z, nsamp = 1e9)
bad("bamscl", Z = sv_z, burn = 2147483647)
ok("bmdsf", D = sv_d, nsamp = 3, burn = 2)
bad("bmdsf", D = matrix(0, 3, 3))
bad("bmdsf", D = sv_d, nsamp = 1e9)
ok("bunfld", D = sv_z, nsamp = 3, burn = 2)
bad("bunfld", D = sv_z, nsamp = 1e9)
ok("cjrirt", votes = sv_v, nsamp = 3, burn = 2, chains = FALSE)
bad("cjrirt", votes = sv_v, seed = -1)
bad("cjrirt", votes = sv_v, nsamp = 1e9)
ok("anomin", votes = cbind(sv_v, 1 - sv_v), nsamp = 2, burn = 1, minvotes = 5)
bad("anomin", votes = sv_v * 2)
bad("anomin", votes = sv_v, ndim = 1e9)
ok("ordirt", Y = sv_z, nsamp = 3, burn = 2)
bad("ordirt", Y = cbind(sv_z, 1))
bad("ordirt", Y = sv_z, nsamp = 1e9)
ok("dynirt", votes = sv_v, period = c(1, 1, 2, 2, 3, 3), nsamp = 3, burn = 2)
bad("dynirt", votes = sv_v, period = 1:2)
bad("dynirt", votes = sv_v, period = c(1, 1, 2, 2, 3, 3), burn = 2147483647)

# Counts over time and across categories
sv_y <- c(52, 48, 61, 55, 47, 50, 66, 70, 58, 49, 45, 53, 57, 51, 63, 60, 49, 52, 71, 74, 61, 50, 47, 55,
          59, 54, 66, 62, 51, 55, 73, 77, 63, 54, 49, 58)
ok("arimfc", y = sv_y, h = 3)
bad("arimfc", y = 1:5)
bad("arimfc", y = sv_y, h = 1e9)
bad("arimfc", y = rep(0.1, 30))
ok("sarfc", y = sv_y, h = 3, order = c(0, 1, 1), seasonal = c(0, 0, 0, 12))
bad("sarfc", y = sv_y, order = c(1, 0, 0), seasonal = c(0, 0, 0, 12))
ok("catcor", units = list(theft = c("a", "b", "b", "c", "d"), fraud = c("b", "c", "c", "d"), arson = c("a", "a", "d")))
bad("catcor", units = list(a = "x"))
ok("ousim", y = sv_y, npath = 2, days = 3)
bad("ousim", y = 1:2)
bad("ousim", y = sv_y, npath = 1, days = 1e7, dt = 1e-9)
ok("fpdens", y = sv_y, ngrid = 8, nstep = 4)
bad("fpdens", y = sv_y, ngrid = 1e9)
ok("lvfit", y = c(410, 455, 430, 470, 520, 498, 515), year = 2017:2023)
bad("lvfit", y = 1:3)
# Exposure and crime: counts by place and day against a daily exposure
# (made-up numbers). The bootstrap cases give a seed: outside R there is
# no generator to draw from without one.
axd <- data.frame(area = rep(letters[1:11], each = 12), day = rep(1:12, 11))
axd$wind <- (axd$day * 47 + seq_len(132) * 13) %% 360
axd$w90 <- as.numeric(axd$wind >= 90 & axd$wind < 180)
axd$w270 <- as.numeric(axd$wind >= 270)
axd$temp <- round(15 + 4 * cos(axd$day / 2) + seq_len(132) %% 5, 2)
axd$pm <- round(20 + 4 * axd$w90 - 3 * axd$w270 + 2 * sin(seq_len(132)) + 0.3 * axd$temp, 2)
axd$pop <- rep(seq(1000, 3000, length.out = 11), each = 12)
axd$crimes <- (axd$day * 3 + seq_len(132)) %% 7 + rep(0:10 %% 3, each = 12) + round(axd$pm / 10)
ok("ppmlfe", data = axd, y = "crimes", x = c("pm", "temp"), fe = c("area", "day"), expo = "pop")
ok("ppmlfe", data = axd, y = "crimes", x = "pm", fe = c("area", "day"), se = "twoway", ssc = "nk")
ok("ppmlfe", data = axd, y = "crimes", x = "pm", fe = "area", se = "nw", time = "day", lag = 2)
ok("ppmlfe", data = axd, y = "crimes", x = c("pm", "temp"), fe = "area", se = "hc1", ssc = "none", level = 0.9)
ok("ppmlfe", data = axd, y = "crimes", x = "pm", fe = c("area", "day"), se = "boot", nboot = 20, seed = 1)
bad("ppmlfe", data = axd, y = "crimes", x = "zz", fe = "area")
bad("ppmlfe", data = axd, y = "crimes", x = "pm", fe = "area", se = "twoway")
ok("ivfe", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"), x = "temp", fe = c("area", "day"))
ok("ivfe", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"), fe = "area", se = "hc1")
ok("ivfe", data = axd, y = "crimes", x = c("pm", "temp"), fe = "area", se = "nw", time = "day")
ok("ivfe", data = axd, y = "crimes", endo = "pm", inst = "w90", fe = "area", se = "boot", nboot = 20, seed = 2)
bad("ivfe", data = axd, y = "crimes", inst = "w90", x = "temp", fe = "area")
ok("ppmlcf", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"), fe = c("area", "day"), nboot = 20,
   seed = 3)
ok("ppmlcf", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"), x = "temp", fe = "area",
   expo = "pop", how = "plug", se = "naive")
bad("ppmlcf", data = axd, y = "crimes", endo = c("pm", "temp"), inst = "w90", fe = "area")
ok("windiv", dir = c(10, 100, 200, 300, NA, 90))
bad("windiv", dir = c(10, 100), width = 360 / 2^27)
ok("windiv", dir = c(10, 100, 200, 300), width = 120, start = 30, area = c("x", "y", "x", "y"), base = 2,
   edges = "paper")
ok("windiv", dir = c(10, 350, 200, 300), sect = c(225, 285))
ok("windiv", dir = c(10, 20, 200, 210, 300, 100), by = c(1, 1, 2, 2, 3, 3), spd = c(1, 2, 1, 1, 3, 1))
ok("windiv", dir = c(10, 20, 200, 210, 300, 100), by = c(1, 1, 2, 2, 3, 3), avg = "mode")
ok("windiv", dir = c(10, 20, 200, 210, 300, 100), by = c(1, 1, 2, 2, 3, 3), share = TRUE,
   area = c("x", "x", "y", "y", "x", "x"))
bad("windiv", dir = c(10, 20), width = 70)
ok("invsn", tlow = c(10, 12, 9, NA), tup = c(12, 11, 9, 8))
ok("invsn", tlow = c(10, 12, 9), tup = c(12, 11, 9.5), sign = "ground", brk = c(-5, 0, 0.5, 5))
ok("invsn", tlow = c(1, 2, 3, 1, 2, 3, 1, NA), tup = c(4, 5, 6, 4, 1, 6, 2, 5), by = c(1, 1, 1, 2, 2, 2, 3, 3))
ok("invsn", tlow = c(1, 2, 3, 1, 2, 3), tup = c(4, 5, 6, 4, 1, 6), by = c("b", "b", "b", "a", "a", "a"),
   rule = "mean", sign = "ground", brk = c(-10, -2, 0, 10))
bad("invsn", tlow = 1, tup = 2, by = c(1, 2))
dws <- data.frame(x = c(0, 50, 0, 20, -40), y = c(0, 0, 0, 10, 30), day = c(1, 1, 2, 2, 3),
                  dir = c(180, 270, 0, 215.5, 95))
dwp <- data.frame(x = c(0, 100, 0, -80, 35), y = c(100, 10, -60, 0, -20))
ok("dwnexp", sx = dws$x, sy = dws$y, sper = dws$day, sdir = dws$dir, rx = dwp$x, ry = dwp$y)
ok("dwnexp", sx = dws$x, sy = dws$y, sper = dws$day, sdir = dws$dir, rx = dwp$x, ry = dwp$y, half = 45, rad = 90,
   from = FALSE, rid = c("a", "b", "c", "d", "e"))
ok("dwnexp", sx = dws$x / 10, sy = dws$y / 10 + 30, sper = dws$day, sdir = dws$dir, rx = dwp$x / 10, ry = dwp$y / 10 + 30,
   rule = "angle", lonlat = TRUE)
ok("dwnexp", sx = dws$x, sy = dws$y, sper = dws$day, sdir = dws$dir, rx = dwp$x, ry = dwp$y, rule = "angle",
   wrap = FALSE, dwt = "div")
ok("dwnexp", sx = c(0, 20), sy = c(0, -10), sper = c(1, 2), sdir = c(200, 120),
   rx = c(-10, 10, 10, -10, 50, 70, 70, 50, -80, 80, 80, -80), ry = c(40, 40, 60, 60, -5, -5, 5, 5, 95, 95, 99, 99),
   rid = rep(1:3, each = 4), rule = "wedge", rad = 100)
bad("dwnexp", sx = 0, sy = 0, sper = NULL, sdir = 0, rx = c(1, 2), ry = c(1, 2), rid = c(1, 1), rule = "wedge")
bad("dwnexp", sx = 0, sy = 0, sper = NULL, sdir = 0, rx = 1, ry = 1, half = 200)
bad("dwnexp", sx = 0, sy = 0, sper = NULL, sdir = 0, rx = 1, ry = 1, rad = 0)
ok("tailwd", dir = c(70, 340, 250, NA, 100), spd = c(4, 4, 4, 2, 3))
ok("tailwd", dir = c(70, 90, 250, 260, NA), spd = c(3, 2, 5, 1, 1), by = c(1, 1, 2, 2, 3), agg = "mean", course = 70)
bad("tailwd", dir = 1, spd = c(1, 2))
bad("invsn", tlow = c(1, 2, 3), tup = c(1, 2))
upd <- data.frame(area = rep(c("w1", "w2", "w3", "e1", "e2", "e3"), 10),
                  ax = rep(c(-100, -300, -600, 150, 350, 700), 10),
                  ay = rep(c(0, 200, -200, 100, -100, 0), 10), day = rep(1:10, each = 6),
                  wind = rep(c(270, 90, 260, 80, 280, 100, 270, 90, 275, 95), each = 6))
upd$y <- (seq_len(60) * 7) %% 11 + (upd$ax > 0) * (upd$wind < 180) * 3
upd$temp <- 10 + upd$day / 2
ups <- data.frame(x0 = 0, y0 = -1000, x1 = 0, y1 = 1000)
ok("updwn", data = upd, y = "y", ax = "ax", ay = "ay", day = "day", wind = "wind", src = ups)
ok("updwn", data = upd, y = "y", ax = "ax", ay = "ay", day = "day", wind = "wind", src = ups, area = "area",
   dist = c(0, 250, 1000), norm = "pois", se = "cluster", clus = "day")
ok("updwn", data = upd, y = "y", ax = "ax", ay = "ay", day = "day", wind = "wind", src = ups, area = "area",
   agg = FALSE, wx = "temp", offs = c(-50, 2000))
bad("updwn", data = upd, y = "y", ax = "ax", ay = "ay", day = "day", wind = "wind", src = data.frame(x0 = 0, y0 = 0))
ok("expjn", expo = data.frame(area = c(1, 1, 2), day = c(1, 2, 1), pm = c(5, 7, 6)),
   outc = data.frame(area = c(1, 2, 2), day = c(1, 1, 3), n = c(2, 4, 1)))
bad("expjn", expo = data.frame(area = 1, day = 1), outc = data.frame(area = 1, day = 1))
skd <- data.frame(place = rep(c("a", "b"), each = 8), day = rep(1:8, 2),
                  plume = c(0, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1),
                  pm = c(8, 20, 7, 9, 25, 11, 8, 10, 6, 7, 18, 8, 9, 30, 7, 22),
                  seas = rep(rep(c("s1", "s2"), each = 4), 2), near = rep(c("b", "a"), each = 8),
                  dens = c("", "light", "", "", "heavy", "light", "", "", "", "", "medium", "", "", "heavy", "", NA))
ok("smkiv", data = skd, plume = "plume", unit = "place", day = "day", pm = "pm", adj = "own", seas = "seas",
   dens = "dens")
ok("smkiv", data = skd, plume = "plume", unit = "place", day = "day", pm = "pm", adj = "nbr", nbr = "near",
   seas = "seas", cf = "median", by = "seas")
ok("smkiv", data = skd, plume = "plume", unit = "place", day = "day")
bad("smkiv", data = skd, plume = "plume", unit = "place", day = "day", pm = "pm")
ok("exprun", data = data.frame(area = rep(1:2, each = 6), day = rep(1:6, 2),
                              smoke = c(0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1)),
   var = "smoke", unit = "area", day = "day")
bad("exprun", data = data.frame(area = 1, day = 1:3, smoke = c(0, 1, 1)), var = "smoke", unit = "area", day = "day",
    cap = 1e9)
ok("exprun", data = data.frame(area = 1, day = c(1, 2, 4, 5), pm = c(30, 40, 50, 10)), var = "pm", unit = "area",
   day = "day", thr = 25, cap = 2)
bad("exprun", data = data.frame(area = 1, day = c(1, 1), s = c(1, 0)), var = "s", unit = "area", day = "day")
# Environmental crime and dynamic panels
en <- data.frame(year = rep(2016:2019, 5), unit = rep(paste0("u", 1:5), each = 4),
                 stringsAsFactors = FALSE)
en$P <- (seq_len(20) * 7) %% 9 + 1
en$S <- (seq_len(20) * 5) %% 11 + 5
en$cases <- (seq_len(20) * 3) %% 13 + 1
ok("lmdix", y = en$cases, q = en[c("P", "S")], unit = en$unit, time = en$year)
ok("lmdix", y = en$cases, q = en[c("P", "S")], unit = en$unit, time = en$year, type = "space",
   ref = "gmean", fnm = c("rate", "dens", "area"))
ok("lmdix", y = replace(en$cases, 3, 0), q = en[c("P", "S")], unit = en$unit, time = en$year,
   at = c(2016, 2018, 2019), zero = 1e-8)
ok("lmdix", y = replace(en$cases, 3, 0), q = cbind(en$P), unit = en$unit, time = en$year, type = "space",
   at = 2017)
bad("lmdix", y = en$cases, q = en[c("P", "S")], unit = en$unit, time = en$year, zero = -1)
bad("lmdix", y = en$cases, q = en[c("P", "S")], unit = en$unit, time = en$year, fnm = "a")
ep <- data.frame(year = rep(2016:2023, 12), prov = rep(sprintf("p%02d", 1:12), each = 8),
                 stringsAsFactors = FALSE)
ep$cper <- 0.02 + 0.01 * sin(seq_len(96))
ep$ur <- 0.02 + 0.005 * cos(seq_len(96) * 1.3)
ep$cases <- 50 + 9000 * ep$cper + 1500 * ep$ur - 4e5 * ep$cper * ep$ur + sin(seq_len(96) * 2.1) +
  as.integer(substring(ep$prov, 2)) + (ep$year - 2016)
ok("fepnl", y = ep$cases, x = ep[c("cper", "ur")], unit = ep$prov, time = ep$year, int = c("cper", "ur"))
ok("fepnl", y = ep$cases, x = ep[c("cper", "ur")], unit = ep$prov, time = ep$year, int = c("cper", "ur"),
   at = c(0.01, 0.02), ssc = "none", level = 0.9)
ok("fepnl", y = ep$cases, x = ep$cper, unit = ep$prov, time = ep$year, clus = rep(1:4, each = 24))
bad("fepnl", y = ep$cases, x = ep[c("cper", "ur")], unit = ep$prov, time = ep$year, int = "cper")
bad("fepnl", y = ep$cases, x = ep[c("cper", "ur")], unit = ep$prov, time = ep$year, at = 1)
ok("pcdtst", x = sin(seq_len(54)) + rep(cos(1:9), 6), unit = rep(paste0("p", 1:6), each = 9),
   time = rep(1:9, 6))
bad("pcdtst", x = rep(1, 4), unit = c(1, 1, 2, 2), time = c(1, 2, 1, 2))
set.seed(41)
gd <- data.frame(t = rep(1:6, 30), u = rep(sprintf("g%02d", 1:30), each = 6), stringsAsFactors = FALSE)
mu <- rnorm(30)
gd$x <- rnorm(180) + mu[match(gd$u, unique(gd$u))]
gd$y <- 0
for (i in seq_len(180)) {
  m <- mu[match(gd$u[i], unique(gd$u))]
  gd$y[i] <- if (gd$t[i] == 1) m / 0.5 + rnorm(1) else 0.5 * gd$y[i - 1] + 0.3 * gd$x[i] + m + rnorm(1)
}
gd <- gd[-c(9, 40), ]
ok("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, lagmax = 3, collapse = TRUE)
ok("gmmpnl", y = gd$y, x = gd["x"], unit = gd$u, time = gd$t, gmm = "x", lagmax = c(3, 2), type = "diff",
   trans = "fod", steps = 2, tdum = FALSE)
ok("gmmpnl", y = gd$y, x = gd["x"], unit = gd$u, time = gd$t, iv = "x", collapse = TRUE, robust = FALSE,
   small = TRUE, hmat = "dpd", ar = 3, level = 0.9)
ok("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, lagmax = 4, steps = 2, robust = FALSE, hmat = "id")
bad("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, steps = 2, wc = TRUE, robust = FALSE)
bad("gmmpnl", y = gd$y, x = gd["x"], unit = gd$u, time = gd$t, gmm = "q")
bad("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, lagmin = 3, lagmax = 2)
ok("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, lagmax = 3, collapse = TRUE, small = TRUE)
ok("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, lagmax = 3, collapse = TRUE, small = TRUE, robust = FALSE)
ok("gmmpnl", y = gd$y, unit = gd$u, time = gd$t, type = "diff", small = TRUE, steps = 2, collapse = TRUE)
bad("gmmpnl", y = c(1, 3, 2, 5, 4, 2), unit = c(1, 1, 1, 1, 1, 2), time = c(1:5, 1), type = "diff",
    tdum = FALSE, small = TRUE)
ok("wsctst", y = gd$y, x = gd["x"], unit = gd$u, time = gd$t)
ok("wsctst", y = gd$y, unit = gd$u, time = gd$t, cons = TRUE, ssc = "none")
bad("wsctst", y = c(1, 3, 2, 4, 5), unit = c(1, 1, 1, 2, 2), time = c(1, 2, 3, 1, 3))
cd <- data.frame(t = rep(1:10, 12), u = rep(sprintf("c%02d", 1:12), each = 10), stringsAsFactors = FALSE)
cf <- sin(cd$t * 1.7)
cg <- (1 + seq_len(12) / 10)[match(cd$u, unique(cd$u))]
cd$x <- cg * cf + cos(seq_len(120) * 0.9)
cd$y <- 1 + 0.7 * cd$x + cg * cf + sin(seq_len(120) * 2.3)
cd$tr <- cd$t
ok("ccepnl", y = cd$y, x = cd["x"], unit = cd$u, time = cd$t)
ok("ccepnl", y = cd$y[-c(5, 33)], x = cd["x"][-c(5, 33), , drop = FALSE], unit = cd$u[-c(5, 33)],
   time = cd$t[-c(5, 33)], w = (1 + match(cd$u, unique(cd$u)) %% 3)[-c(5, 33)], level = 0.9)
ok("ccepnl", y = cd$y, x = cd["x"], unit = cd$u, time = cd$t, fac = cd["tr"], cons = FALSE, avg = "x")
bad("ccepnl", y = cd$y, x = ave(cd$x, cd$t), unit = cd$u, time = cd$t)
bad("ccepnl", y = cd$y, x = cd["x"], unit = cd$u, time = cd$t, avg = "q")
ok("ccepnl", y = cd$y, x = cd["x"], unit = cd$u, time = cd$t, pvar = "bsw")
ok("ccepnl", y = cd$y, x = cd["x"], unit = cd$u, time = cd$t, pool = FALSE)
dd <- data.frame(t = rep(1:24, 10), u = rep(sprintf("d%02d", 1:10), each = 24), stringsAsFactors = FALSE)
df <- sin(dd$t * 1.3)
dg <- (1 + seq_len(10) / 7)[match(dd$u, unique(dd$u))]
dd$x <- dg * df + cos(seq_len(240)^2 * 0.37)
dd$y <- 0
for (r in which(dd$t > 1)) {
  dd$y[r] <- 0.4 * dd$y[r - 1] + 0.7 * dd$x[r] + dg[r] * df[r] + sin(r^2 * 0.11)
}
ok("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, ylag = 1)
ok("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, ylag = 1, lags = 1, bias = "jack")
ok("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, ylag = 1, lags = 1, bias = "rma")
bad("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, ylag = 1, pool = TRUE)
bad("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, bias = "jack", pool = TRUE)
bad("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, ylag = 1, lags = 23)
bad("ccepnl", y = dd$y, x = dd["x"], unit = dd$u, time = dd$t, pvar = "bsw", w = (seq_along(dd$y) - 1) %/% 24 + 1)
# Regression discontinuity
set.seed(20261120)
rdd <- data.frame(r = round(runif(300, -1, 1), 4))
rdd$y <- round(1 + 0.8 * rdd$r + 2 * rdd$r^2 * (rdd$r >= 0) +
                 0.5 * (rdd$r >= 0) + rnorm(300, 0, 0.3), 4)
rdd$t <- as.numeric(runif(300) < 0.2 + 0.6 * (rdd$r >= 0))
rdd$g <- sample(1:30, 300, TRUE)
ok("rdloc", data = rdd, y = "y", run = "r")
ok("rdloc", data = rdd, y = "y", run = "r", trt = "t", bw = "user", h = 0.6,
   kern = "epa", ci = "conventional")
ok("rdloc", data = rdd, y = "y", run = "r", bw = "user", h = 0.5,
   kern = "uni", se = "cluster", clus = "g", ssc = "none",
   plac = c(-0.5, 0.5), level = 0.9)
bad("rdloc", data = rdd, y = "y", run = "r", bw = "user")
bad("rdloc", data = rdd, y = "y", run = "r", bw = "user", h = 0.001)
bad("rdloc", data = rdd, y = "y", run = "r", plac = 0)
bad("rdloc", data = rdd, y = "y", run = "r", cut = 5)
ok("rddens", r = rdd$r, b = 0.05, h = 0.5)
ok("rddens", r = rdd$r, cut = 0.1, level = 0.9)
bad("rddens", r = c(1, NA))
# The exponential-mean GMM and the functions that take an estimator
# function. Outside R `est` names the estimator: these calls give R the
# function of that name.
for (f in c("runfx", "hetfx", "ccross", "jntx", "dosebin", "typefx",
            "plcday")) {
  assign(f, local({
    g <- getExportedValue("crimlab", f)
    function(...) {
      a <- list(...)
      if (is.character(a$est)) a$est <- match.fun(a$est)
      do.call(g, a)
    }
  }))
}
ok("expgmm", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"),
   x = "temp", fe = "area", se = "hc1")
ok("expgmm", data = axd, y = "crimes", endo = "pm", inst = c("w90", "w270"),
   fe = c("area", "day"), se = "hc1", step = "one", err = "add", ssc = "none",
   level = 0.9)
ok("expgmm", data = axd, y = "crimes", endo = "pm", inst = "w90",
   expo = "pop", se = "hc1", step = "iter")
bad("expgmm", data = axd, y = "crimes", endo = "pm", inst = "pm")
bad("expgmm", data = axd, y = "crimes", endo = "pm", inst = "w90")
ok("runfx", data = axd, y = "crimes", var = "pm", unit = "area", day = "day",
   thr = 22, fe = c("area", "day"))
ok("runfx", data = axd, y = "crimes", var = "pm", unit = "area", day = "day",
   thr = 21, cap = 2, ctrl = "temp", fe = "area", se = "hc1")
bad("runfx", data = axd, y = "crimes", var = "pm", unit = "area", day = "day",
    est = "sum", fe = "area")
axd$g <- rep(c("n", "s", "s"), 44)
ok("hetfx", data = axd, var = "pm", mod = "g", y = "crimes",
   x = c("pm", "temp"), fe = "area")
ok("hetfx", data = axd, var = "pm", mod = "g", how = "split", y = "crimes",
   x = "pm", fe = "area", se = "hc1")
ok("hetfx", data = axd, var = "pm", mod = "g", est = "ivfe", y = "crimes",
   endo = "pm", inst = "w90", x = "temp", fe = "area", se = "hc1")
bad("hetfx", data = axd, var = "pm", mod = "area", y = "crimes", x = "pm",
    fe = "day", how = "x")
# eight weeks for four areas: a stratum of year, month and weekday then
# holds several days, which twelve days running do not give
axc <- data.frame(area = rep(letters[1:4], each = 56), day = rep(1:56, 4))
axc$date <- as.Date("2024-01-07") + axc$day
ic <- seq_len(nrow(axc))
axc$temp <- round(15 + 4 * cos(axc$day / 2) + ic %% 5, 2)
axc$pm <- round(20 + 3 * sin(ic) + 2 * cos(ic / 3) + 0.3 * axc$temp, 2)
axc$crimes <- (axc$day * 3 + ic) %% 7 + rep(0:3 %% 3, each = 56) +
  round(axc$pm / 10)
ok("ccross", data = axc, y = "crimes", unit = "area", day = "date", x = "pm",
   se = "hc1")
ok("ccross", data = axc, y = "crimes", unit = "area", day = "date",
   strat = "ym", ma = "pm", k = 3, x = "pm", ns = c(temp = 3), se = "hc1")
ok("ccross", data = axc, y = "crimes", unit = "area", day = "date",
   lagv = "temp", lag = 1, x = c("pm", "temp"), k = 1)
bad("ccross", data = axc, y = "crimes", unit = "area", day = "date", k = 0)
ok("jntx", data = axd, y = "crimes", x1 = "pm", x2 = "temp", fe = "area",
   se = "hc1")
ok("jntx", data = axd, y = "crimes", x1 = "pm", x2 = "temp", form = "linear",
   ref = c(20, 15), ngrid = 3, p2 = c(0.25, 0.75), fe = "area")
bad("jntx", data = axd, y = "crimes", x1 = "pm", x2 = "temp", p1 = 0.5,
    fe = "area")
ok("dosebin", data = axd, y = "crimes", var = "pm", fe = "area")
ok("dosebin", data = axd, y = "crimes", var = "pm", brk = c(0, 22, 25, 40),
   ref = 2, ctrl = "temp", fe = "area", se = "hc1")
bad("dosebin", data = axd, y = "crimes", var = "pm", brk = c(1, 2),
    fe = "area")
axd$c2 <- (axd$day * 5 + seq_len(132)) %% 4 + round(axd$temp / 6)
ok("typefx", data = axd, ys = c("crimes", "c2"), term = "pm",
   x = c("pm", "temp"), fe = "area")
ok("typefx", data = axd, ys = c("crimes", "c2"), term = "pm", adj = "BH",
   x = "pm", fe = c("area", "day"), se = "hc1")
bad("typefx", data = axd, ys = "crimes", term = c("a", "b"), x = "pm",
    fe = "area")
bad("typefx", data = axd, ys = "crimes", term = "pm", est = "ivfe",
    x = "pm", fe = "area", maxit = 5)
ok("plcday", data = axd, var = "pm", unit = "area", day = "day",
   shifts = c(-1, 1), y = "crimes", x = c("pm", "temp"), fe = "area")
ok("plcday", data = axd, var = "pm", unit = "area", day = "day", shifts = 2,
   est = "ivfe", with = "w90", y = "crimes", endo = "pm", inst = "w90",
   x = "temp", fe = "area", se = "hc1")
bad("plcday", data = axd, var = "pm", unit = "area", day = "day", shifts = 0,
    y = "crimes", x = "pm", fe = "area")
# The palettes, and the rows the plots draw (drawn on no device here)
grDevices::pdf(NULL)
ok("blpal")
ok("blpal", name = "sequential", n = 7)
ok("blpal", name = "grey", n = 7)
ok("blpal", name = "diverging", n = 1)
bad("blpal", name = "x")
bad("blpal", n = 0)
ok("yoypals")
ok("plotcdfs", x = c(3, 1, 2, NA, 5, 4), grp = c("b", "a", "b", "a", "b", "a"))
bad("plotcdfs", x = "a")
ok("plotlrnz", x = c(5, 1, 3, 0, 8))
ok("plotlrnz", x = data.frame(pop = c(1, 0.5), share = c(1, 0.2)))
bad("plotlrnz", x = c(-1, 2))
ok("plotgain", x = c(5, 1, 3, 0, 8, 2), grp = rep(c("u", "v"), 3))
bad("plotgain", x = c(5, -1, 3))
ok("plotmult", y = c(3, 4, 5, 6), prd = c(2, 1, 2, 1), grp = c("a", "a", "b", "b"))
ok("plotmult", y = c(3, 4, NA), prd = c("q2", "q1", "q3"))
ok("plothrz", y = c(3, 4, 5, 6, 9, 1), prd = rep(1:3, 2),
   grp = rep(c("a", "b"), each = 3), bands = 2)
bad("plothrz", y = 1:3, prd = 1:3, bands = 9)
ok("plotbump", y = c(3, 4, 5, 5), prd = c(1, 1, 2, 2), grp = c("a", "b", "a", "b"))
ok("plotbump", y = c(3, 4, 5, 6), prd = c(1, 1, 2, 2), grp = c("a", "b", "a", "b"),
   desc = FALSE)
bad("plotbump", y = 1:2, prd = 1:2, grp = c("a", "b"), desc = "x")
ok("plotfnl", obs = c(5, 20, 3, 12), expd = c(6, 10, 4, 11),
   area = c("p", "q", "r", "s"))
ok("plotfnl", obs = c(5, 20, NA, 12), expd = c(6, 10, 4, 11), lvls = 0.9,
   target = 1.2)
bad("plotfnl", obs = 5, expd = 6, lvls = 1)
ok("plotcal", x = as.Date("2024-02-27") + 0:5, y = c(1, 4, 2, 8, 5, 7))
bad("plotcal", x = as.Date("2024-02-27"), y = "a")
invisible(grDevices::dev.off())
out <- file.path("cli", "test.tsv")
old <- readLines(out, encoding = "UTF-8")
mark <- grep("^# made by dev/clref.R", old)
keep <- if (length(mark)) old[seq_len(mark[1L] - 1L)] else old
writeLines(c(keep, "# made by dev/clref.R: do not edit the lines below", cases), out, useBytes = TRUE)
cat(length(cases), "cases written to", out, "\n")
