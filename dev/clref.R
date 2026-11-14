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
  paste0('"', gsub('"', '\\\\"', gsub("\\\\", "\\\\\\\\", enc2utf8(s))), '"')
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
  paste0("{", paste(sprintf("%s:%s", esc(names(l)), vapply(l, js, "")), collapse = ","), "}")
}
plain <- function(x) {
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
ok <- function(.f, ...) {
  args <- list(...)
  add(.f, args, js(do.call(.f, args)))
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
set.seed(9)
tby <- round(rnorm(30, 10), 3)
tbx <- data.frame(temp = round(runif(30, 5, 30), 2))
ok("tpbhm", y = tby, x = tbx)
ok("tpbhm", y = tby, x = tbx, grp = rep(c("p", "q", "r"), 10), pm = c(0, 0.1), psd = 5, level = 0.9)
bad("tpbhm", y = tby, x = tbx[1:5, , drop = FALSE])

# Learners
set.seed(10)
lx <- matrix(round(rnorm(60), 3), 20, 3)
ly <- round(lx[, 1] - 0.5 * lx[, 2] + rnorm(20, sd = 0.3), 3)
lb <- rbinom(20, 1, plogis(0.8 * lx[, 1]))
ok("enpath", x = lx, y = ly, nlam = 5)
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
sy2 <- c(3, 4, 3, 5, 4, 9, 10, 8, 11, 9)
ok("stpcp", y = sy2, nperm = 99, seed = 1)
ok("stpcp", y = sy2, x = 2001:2010, minseg = 3, nperm = 49, seed = 2)
bad("stpcp", y = sy2, x = 1:3, seed = 1)

out <- file.path("cli", "test.tsv")
old <- readLines(out, encoding = "UTF-8")
mark <- grep("^# made by dev/clref.R", old)
keep <- if (length(mark)) old[seq_len(mark[1L] - 1L)] else old
writeLines(c(keep, "# made by dev/clref.R: do not edit the lines below", cases), out, useBytes = TRUE)
cat(length(cases), "cases written to", out, "\n")
