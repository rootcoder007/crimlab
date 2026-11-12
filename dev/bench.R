# Speed of every exported function against every package found that
# computes the same quantity, on fixed seeded inputs. A pair is timed only
# after the two results agree to the stated relative tolerance; a pair
# that disagrees is reported as such and fails the run. Where no package
# computes the quantity, the counterpart is a plain base-R implementation
# of the same formula ("base R"). Calls are timed as a user makes them,
# argument checks included; graphics go to a null device.
#
#   Rscript dev/bench.R [OUTDIR] [fn ...]
#
# writes OUTDIR/bench.csv and OUTDIR/bench.md (OUTDIR defaults to
# inst/bench) and, when every function was run, the list of functions in
# README.md between the bench markers. Exit status 1 on a disagreement or
# an error, never on a slow result. Columns: fn, case (which arguments,
# where one function is timed two ways), ref (package::function),
# n and of (what n counts), crimlab and other (median seconds per call),
# ratio (other / crimlab: above 1, crimlab is faster), agree, tol and scale
# (relative to each value or to the largest), refver (the counterpart's
# version), r, cpu, date.
#
# Counterparts: rmorie (r-universe, not CRAN) for the disparity, policing
# and corrections functions; Familias and forrel for kinidx, patidx and
# (Familias) gtmp; forensim for mixlr; comparison for cllr, cllrmin, lrcal
# and eceplt. forensim declares tkrplot, which needs Tk headers; only
# forensim::LR is used.
a0 <- commandArgs(TRUE)
out <- if (length(a0)) a0[1] else "inst/bench"
only <- a0[-1]
suppressMessages(library(crimlab))
grDevices::pdf(NULL)
have <- function(p) requireNamespace(p, quietly = TRUE)
nm <- function(k, x) stats::setNames(x, paste0("g", seq_len(k)))

cases <- list()
# va and vb pick the values compared from the two results (and may read
# the inputs, as a second argument); tol is relative to the larger of the
# two values, or with sc to the largest value of either result (for terms
# that cancel)
add <- function(fn, ref, n, of, mk, a, b, va, vb, tol = 1e-9, sc = FALSE,
                case = "") {
  cases[[length(cases) + 1L]] <<- list(
    fn = fn, case = case, ref = ref, n = n, of = of, mk = mk, a = substitute(a),
    b = substitute(b), va = va, vb = vb, tol = tol, sc = sc
  )
}
pick <- function(...) {
  k <- c(...)
  function(r) unlist(lapply(k, function(j) r[[j]]))
}

## disparities ------------------------------------------------------------

add("vodtst", "rmorie::morie_veil_of_darkness", c(5, 1e3, 1e5), "groups",
  function(n) list(d = nm(n, stats::runif(n, 10, 1e3)),
                   e = nm(n, stats::runif(n, 10, 1e3))),
  vodtst(d, e), rmorie::morie_veil_of_darkness(d, e),
  pick("share_day", "share_night", "ratio"),
  pick("share_day", "share_night", "bias_ratio"))

add("orrisk", "rmorie::morie_relative_risk_from_or", 1, "calls",
  function(n) list(),
  orrisk(3, 0.4, 0.5), rmorie::morie_relative_risk_from_or(3, 0.4, 0.5),
  function(r) r[c("exposed", "unexposed", "rr")],
  function(r) r[["risks"]],
  tol = 1e-6) # rmorie finds the root numerically

add("outtst", "rmorie::morie_outcome_test", c(1e3, 1e5, 1e6), "people",
  function(n) list(q = stats::runif(n), w = stats::runif(n)),
  outtst(q, 0.5, w), rmorie::morie_outcome_test(list(a = list(
    q = q, w = w, t = 0.5
  ))),
  pick("stop", "hit"), pick("stop_rate", "hit_rate"))

add("dspbm", "rmorie::morie_disparity_benchmark", c(5, 1e3, 1e5), "groups",
  function(n) list(p = nm(n, stats::runif(n, 1e3, 1e5)),
                   cn = nm(n, stats::runif(n, 10, 1e3)),
                   f = nm(n, stats::runif(n, 1, 50))),
  dspbm(p, cn, f), rmorie::morie_disparity_benchmark(p, cn, f, "g1"),
  pick("cd", "fd", "rd"),
  pick("contact_disparity", "force_given_contact_disparity",
      "resident_disparity"))

add("dyadrt", "rmorie::morie_interracial_rates", c(5, 50, 300), "groups",
  function(n) {
    p <- nm(n, stats::runif(n, 1e3, 1e5))
    m <- matrix(stats::runif(n * n, 0, 100), n)
    v <- c(m)
    names(v) <- paste0(names(p)[row(m)], "_on_", names(p)[base::col(m)])
    list(m = m, p = p, v = v)
  },
  dyadrt(m, p), rmorie::morie_interracial_rates(v, p),
  pick("rate_off", "rate_pair", "null_off", "ratio"),
  pick("rate_per_offender_group", "rate_per_pair_exposure",
      "null_rate_per_offender_group", "ratio_to_null"))

## forensics --------------------------------------------------------------

# six alleles a locus; g1, g2 and the trio drawn from them
dna <- function(n) {
  fq <- t(replicate(n, {
    f <- stats::runif(6)
    f / sum(f)
  }))
  s <- function() t(replicate(n, sample(6, 2, TRUE)))
  g1 <- s()
  g2 <- s()
  mo <- s()
  fa <- s()
  ch <- cbind(mo[cbind(seq_len(n), sample(2, n, TRUE))],
              fa[cbind(seq_len(n), sample(2, n, TRUE))])
  fr <- function(g) matrix(fq[cbind(rep(seq_len(n), 2), c(g))], n)
  list(fq = fq, g1 = g1, g2 = g2, f1 = fr(g1), mo = mo, fa = fa, ch = ch,
       fc = fr(ch))
}
# Familias: likelihoods per locus of pedigrees for persons ids with
# genotypes gs (a list of n x 2 matrices)
fam <- function(fq, peds, ids, gs) {
  loc <- lapply(seq_len(nrow(fq)), function(i) {
    Familias::FamiliasLocus(fq[i, ], allelenames = 1:6, name = paste0("L", i))
  })
  d <- do.call(cbind, lapply(seq_len(nrow(fq)), function(i) {
    t(vapply(gs, function(g) g[i, ], numeric(2)))
  }))
  d <- as.data.frame(d, row.names = ids)
  names(d) <- paste0("L", rep(seq_len(nrow(fq)), each = 2), c(".1", ".2"))
  r <- Familias::FamiliasPosterior(peds, loc, d, ref = length(peds))
  r[["likelihoodsPerSystem"]]
}
fped <- function(id, dad, mom, sex) {
  Familias::FamiliasPedigree(id = id, dadid = dad, momid = mom, sex = sex)
}
famki <- function(fq, g1, g2) {
  fs <- fped(c("a", "b", "f", "m"), c("f", "f", NA, NA), c("m", "m", NA, NA),
             c("male", "male", "male", "female"))
  un <- fped(c("a", "b"), c(NA, NA), c(NA, NA), c("male", "male"))
  l <- fam(fq, list(fs, un), c("a", "b"), list(g1, g2))
  sum(log10(l[, 1] / l[, 2]))
}
fampi <- function(fq, ch, fa, mo) {
  ids <- c("mo", "ch", "af")
  pf <- fped(ids, c(NA, "af", NA), c(NA, "mo", NA),
             c("female", "female", "male"))
  pu <- fped(c(ids, "tf"), c(NA, "tf", NA, NA), c(NA, "mo", NA, NA),
             c("female", "female", "male", "male"))
  l <- fam(fq, list(pf, pu), ids, list(mo, ch, fa))
  sum(log10(l[, 1] / l[, 2]))
}
# forrel: genotypes set on the first pedigree, moved to the alternatives
frl <- function(fq, h1, h2, gs) {
  am <- matrix(0, length(labels(h1)), 2 * nrow(fq),
               dimnames = list(labels(h1), NULL))
  for (id in names(gs)) am[id, ] <- c(t(gs[[id]]))
  la <- lapply(seq_len(nrow(fq)), function(i) {
    list(alleles = 1:6, afreq = fq[i, ], name = paste0("L", i))
  })
  h1 <- pedtools::setMarkers(h1, alleleMatrix = am, locusAttributes = la)
  h2 <- pedtools::transferMarkers(h1, h2)
  r <- forrel::kinshipLR(h1, h2)
  r[["lnLRtotal"]][[1L]] / log(10)
}
frlki <- function(fq, g1, g2) {
  frl(fq, pedtools::nuclearPed(2),
      list(pedtools::singleton(3), pedtools::singleton(4)),
      list("3" = g1, "4" = g2))
}
frlpi <- function(fq, ch, fa, mo) {
  frl(fq, pedtools::nuclearPed(father = "af", mother = "mo", children = "ch"),
      list(pedtools::nuclearPed(father = "tf", mother = "mo",
                                children = "ch"),
           pedtools::singleton("af")),
      list(ch = ch, af = fa, mo = mo))
}
famgt <- function(p, q, hom, th) {
  n <- length(p)
  fq <- cbind(p, q, 1 - p - q)
  loc <- lapply(seq_len(n), function(i) {
    Familias::FamiliasLocus(fq[i, ], allelenames = c("A", "B", "C"),
                            name = paste0("L", i))
  })
  g <- ifelse(hom, "A", "B")
  un <- fped(c("a", "b"), c(NA, NA), c(NA, NA), c("male", "male"))
  one <- fped("a", NA, NA, "male")
  dd <- function(k) {
    d <- as.data.frame(matrix(rep(c(rbind("A", g)), each = k), k),
                       row.names = c("a", "b")[seq_len(k)])
    names(d) <- paste0("L", rep(seq_len(n), each = 2), c(".1", ".2"))
    d
  }
  l2 <- Familias::FamiliasPosterior(list(un), loc, dd(2), ref = 1,
                                    kinship = th)[["likelihoodsPerSystem"]]
  l1 <- Familias::FamiliasPosterior(list(one), loc, dd(1), ref = 1,
                                    kinship = th)[["likelihoodsPerSystem"]]
  sum(log10(l2 / l1))
}
gmk <- function(n) {
  list(p = stats::runif(n, 0.02, 0.3), q = stats::runif(n, 0.02, 0.3),
       h = stats::runif(n) < 0.2)
}
add("gtmp", "Familias::FamiliasPosterior", c(20, 200), "loci", gmk,
  gtmp(p, q, h, 0.01), famgt(p, q, h, 0.01),
  function(r) r[["log10mp"]], identity)

add("kinidx", "Familias::FamiliasPosterior", c(20, 200), "loci", dna,
  kinidx(g1, g2, f1, "sibling"), famki(fq, g1, g2),
  function(r) r[["log10idx"]], identity)
add("kinidx", "forrel::kinshipLR", c(20, 200), "loci", dna,
  kinidx(g1, g2, f1, "sibling"), frlki(fq, g1, g2),
  function(r) r[["log10idx"]], identity)
add("patidx", "Familias::FamiliasPosterior", c(20, 200), "loci", dna,
  patidx(ch, fa, fc, mo), fampi(fq, ch, fa, mo),
  function(r) r[["log10idx"]], identity)
add("patidx", "forrel::kinshipLR", c(20, 200), "loci", dna,
  patidx(ch, fa, fc, mo), frlpi(fq, ch, fa, mo),
  function(r) r[["log10idx"]], identity)

# a mixture of n alleles: the victim explains alleles 1 and 2; under the
# first account the suspect explains the rest, under the second
# (n - 2) / 2 unknowns do
add("mixlr", "forensim::LR", c(4, 6, 8), "alleles",
  function(n) {
    p <- stats::setNames(stats::runif(n, 0.02, 0.08), seq_len(n))
    list(p = p, al = names(p), k = (n - 2) / 2,
         fr = c(p, "99" = 1 - sum(p)))
  },
  mixlr(p, character(0), 0, al[-(1:2)], k),
  forensim::LR(Repliste = seq_along(p), Tp = seq_along(p), Td = 1:2,
               Vp = 0, Vd = 0, xp = 0, xd = k, theta = 0,
               prDHet = rep(0, 8), prDHom = rep(0, 8), prC = 0, freq = fr),
  pick("num", "den", "lr"), pick("num", "deno", "LR"))

lrs <- function(n) {
  list(s = exp(stats::rnorm(n / 2, 1, 1.5)), d = exp(stats::rnorm(n / 2, -1, 1.5)))
}
add("cllr", "comparison::calc.ece", c(1e3, 1e4), "ratios", lrs,
  cllr(s, d), comparison::calc.ece(s, d, 0.5),
  identity, function(r) r[["ece"]])
add("cllr", "base R", c(1e3, 1e5, 1e6), "ratios", lrs,
  cllr(s, d), (mean(log2(1 + 1 / s)) + mean(log2(1 + d))) / 2,
  identity, identity)
add("cllrmin", "comparison::calc.ece", c(1e3, 1e4), "ratios", lrs,
  cllrmin(s, d), comparison::calc.ece(s, d, 0.5),
  pick("cllr", "min"), pick("ece", "ece.cal"))
# comparison pools with isotone::gpava, whose time grows with the square of
# n (about 6 minutes at 1e5 on the reference machine): 1e3 and 1e4 only.
# calibrate.set returns the ratios sorted, not in the order given
add("lrcal", "comparison::calibrate.set", c(1e3, 1e4), "ratios", lrs,
  lrcal(s, d), comparison::calibrate.set(s, d),
  function(r, e) c(r[["same"]][order(e$s)], r[["diff"]][order(e$d)]),
  pick("LR.cal.ss", "LR.cal.ds"))
add("eceplt", "comparison::calc.ece", c(1e3, 1e4), "ratios", lrs,
  eceplt(s, d), {
    e <- comparison::calc.ece(s, d, 1 / (1 + 10^-seq(-2.5, 2.5, by = 0.05)))
    graphics::plot(e)
    e
  },
  pick("ece", "min", "ref"), pick("ece", "ece.cal", "ece.null"))

add("lrpost", "base R", c(1e3, 1e5, 1e6), "ratios",
  function(n) list(l = exp(stats::rnorm(n, 0, 5))),
  lrpost(l, 1e-3), {
    o <- l * 1e-3 / (1 - 1e-3)
    o / (1 + o)
  }, identity, identity)
add("dbsrch", "base R", 1, "calls", function(n) list(),
  dbsrch(1e-6, 1e5, 3e7),
  c(np = 1e-6 * 1e5, any = -expm1(1e5 * log1p(-1e-6)),
    post = 1 / (1 + (3e7 - 1e5) * 1e-6)),
  identity, identity)
add("tipplt", "base R", c(1e3, 1e5, 1e6), "ratios", lrs, tipplt(s, d), {
  f <- function(x) {
    x <- sort(log10(x))
    data.frame(log10lr = x, share = rev(seq_along(x)) / length(x))
  }
  r <- list(same = f(s), diff = f(d))
  graphics::plot(range(r$same$log10lr, r$diff$log10lr), c(0, 1), type = "n",
                 xlab = "log10 likelihood ratio",
                 ylab = "share of ratios at or above")
  graphics::lines(r$same$log10lr, r$same$share, type = "s", lty = 1)
  graphics::lines(r$diff$log10lr, r$diff$share, type = "s", lty = 2)
  graphics::abline(v = 0, lty = 3)
  graphics::legend("bottomleft", c("same source", "different source"),
                   lty = 1:2, bty = "n")
  r
}, function(r) unlist(r), function(r) unlist(r))

## policing ---------------------------------------------------------------

# Kitagawa's split, the mix term of each part centred on the mean of the
# two aggregate rates. A part's mix term is a difference of nearly equal
# rates; with many parts it can be 1e-6 of the largest term, and the two
# orders of operation then differ in its last digits: tolerance against
# the largest term.
add("clrsmp", "base R", c(5, 1e3, 1e5), "offence types",
  function(n) {
    ra <- stats::runif(n, 10, 1e3)
    rb <- stats::runif(n, 10, 1e3)
    list(ra = ra, rb = rb, ca = ra * stats::runif(n), cb = rb * stats::runif(n))
  },
  clrsmp(ca, ra, cb, rb), {
    a <- ca / ra
    b <- cb / rb
    sa <- ra / sum(ra)
    sb <- rb / sum(rb)
    ta <- sum(ca) / sum(ra)
    tb <- sum(cb) / sum(rb)
    re <- (b - a) * (sa + sb) / 2
    me <- (sb - sa) * ((a + b) / 2 - (ta + tb) / 2)
    c(re, me, ta, tb, tb - ta, sum(re), sum(me))
  },
  function(r) c(r$parts$rateeff, r$parts$mixeff, r$total[1:5]), identity,
  sc = TRUE)

add("expbnd", "rmorie::morie_disparity_exposure_bounds", c(5, 1e3, 1e5),
  "groups",
  function(n) list(y = nm(n, stats::runif(n, 10, 1e3)),
                   m = nm(n, stats::runif(n, 1e3, 1e4))),
  expbnd(y, m, 1.5), rmorie::morie_disparity_exposure_bounds(y, m, 1.5, "g1"),
  pick("rate", "lo", "hi", "rlo", "rhi"),
  pick("proxy_rate", "rate_lower", "rate_upper", "ratio_lower",
      "ratio_upper"))

add("misrd", "rmorie::morie_misclassified_rd", 1, "calls", function(n) list(),
  misrd(0.3, 0.1, 0.4, 0.8, 0.9),
  rmorie::morie_misclassified_rd(0.3, 0.1, 0.4, 0.8, 0.9),
  pick("true", "measured", "factor", "d1", "d0"),
  pick("true_rd", "measured_rd", "factor", "d1", "d0"))

add("rtmhs", "rmorie::morie_regression_to_mean", c(1e3, 1e5, 1e6), "places",
  function(n) {
    mu <- stats::rgamma(n, 2, 0.5)
    x1 <- stats::rpois(n, mu)
    list(x1 = x1, x2 = stats::rpois(n, mu),
         cut = stats::quantile(x1, 0.9, names = FALSE))
  },
  rtmhs(x1, x2, cut), rmorie::morie_regression_to_mean(x1, x2, cut),
  pick("nsel", "change", "mirror", "sym", "excess", "low", "lowsym"),
  pick("n_selected", "selected_change", "mirror_change", "symmetrised_change",
      "excess_over_symmetry", "low_selected_change",
      "low_symmetrised_change"))

## corrections ------------------------------------------------------------

add("incap", "rmorie::morie_incapacitation", c(5, 1e3, 1e5), "groups",
  function(n) {
    sh <- stats::runif(n)
    list(l = stats::runif(n, 0.5, 20), sh = sh / sum(sh))
  },
  incap(l, 0.1, 2, share = sh),
  rmorie::morie_incapacitation(l, 0.1, 2, shares = sh), case = "share",
  function(r) {
    k <- seq_len(nrow(r) - 1L)
    c(r$rate[k], r$prev[k], r$marg[k], r$prev[nrow(r)])
  },
  function(r) {
    c(r$rate, r$prevented_share, r$marginal_prevention,
      attr(r, "aggregate")[["prevented_share"]])
  })

add("incap", "rmorie::morie_incapacitation", c(5, 1e3, 1e5), "groups",
  function(n) list(l = stats::runif(n, 0.5, 20)),
  incap(l, 0.1, 2), rmorie::morie_incapacitation(l, 0.1, 2),
  pick("rate", "prev", "marg"),
  pick("rate", "prevented_share", "marginal_prevention"), case = "no share")

add("incapr", "rmorie::morie_incapacitation_career", c(10, 1e3, 1e5),
  "career periods",
  function(n) list(l = sort(stats::runif(n, 0, 20), decreasing = TRUE),
                   s = floor(n / 2)),
  incapr(l, 2, s, 0.25), rmorie::morie_incapacitation_career(l, 2, s, 0.25),
  pick("prev", "net", "upper", "lower", "later"),
  pick("prevented", "prevented_net", "upper", "lower", "later"))

add("sellab", "rmorie::morie_selective_labels", c(1e3, 1e5, 1e6), "defendants",
  function(n) {
    k <- stats::runif(n)
    list(y = stats::rbinom(n, 1, k), rel = k < 0.7, rule = k < 0.8)
  },
  sellab(y, rel, rule), rmorie::morie_selective_labels(y, rel, rule),
  function(r) r[c("obs", "unobs", "lower", "upper", "naive")],
  function(r) c(r$observed_rate, r$rule_share_unobserved, r$bounds,
                r$naive_rate))

add("losbnd", "rmorie::morie_backlog_censoring", c(1e3, 1e5, 1e6), "stays",
  function(n) list(d = stats::rexp(n / 2, 1 / 60), o = stats::rexp(n / 2, 1 / 90)),
  losbnd(d, o), rmorie::morie_backlog_censoring(d, o),
  pick("mean", "age", "lower", "bias"),
  pick("disposed_mean", "pending_age", "lower_bound", "bias_lower"))

add("sntbnd", "rmorie::morie_sentence_effect_bounds", c(1e3, 1e5, 1e6),
  "people",
  function(n) {
    z <- sample(c("community", "custody"), n, TRUE, c(0.6, 0.4))
    list(z = z, y = stats::rbinom(n, 1, ifelse(z == "custody", 0.55, 0.35)))
  },
  sntbnd(y, z), rmorie::morie_sentence_effect_bounds(y, z),
  function(r) c(r$lower, r$upper, r$observed[3]),
  function(r) {
    c(r$outcome_bounds[, "lower"], r$ate_bounds[["lower"]],
      r$outcome_bounds[, "upper"], r$ate_bounds[["upper"]],
      r$naive_difference)
  })

add("hzdsel", "rmorie::morie_hazard_selection", 1, "calls", function(n) list(),
  hzdsel(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9)),
  rmorie::morie_hazard_selection(0.5, 0.5, 0.1, c(0.5, 0.8), c(0.9, 0.9)),
  pick("share0", "share1", "haz0", "haz1", "ratio"),
  pick("surviving_high_share", "period2_hazard", "period2_hazard_ratio"))

add("proof", "base R", 1, "calls", function(n) list(tab = proof()),
  proof(c("cllr", "gtmp")), tab[tab$fn %in% c("cllr", "gtmp"), ],
  function(r) unlist(r), function(r) unlist(r))

## run --------------------------------------------------------------------

same <- function(x, y, tol, sc) {
  x <- unname(unlist(x))
  y <- unname(unlist(y))
  if (!is.numeric(x) || !is.numeric(y)) return(identical(x, y))
  k <- !is.na(x)
  if (length(x) != length(y) || !identical(k, !is.na(y))) return(FALSE)
  m <- if (sc) max(abs(c(x[k], y[k]))) else pmax(abs(x[k]), abs(y[k]))
  isTRUE(all(x[k] == y[k] | abs(x[k] - y[k]) <= tol * m))
}
val <- function(f, r, e) if (length(formals(f)) == 2L) f(r, e) else f(r)
cpu <- function() {
  m <- tryCatch(system2("lscpu", stdout = TRUE, stderr = FALSE),
                error = function(e) character(0), warning = function(w) character(0))
  m <- sub("^Model name: *", "", grep("^Model name:", m, value = TRUE))
  if (!length(m)) {
    m <- tryCatch(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                          stdout = TRUE, stderr = FALSE),
                  error = function(e) character(0),
                  warning = function(w) character(0))
  }
  if (length(m)) paste(m[1], Sys.info()[["machine"]]) else
    Sys.info()[["machine"]]
}
env <- c(r = paste(R.version$major, R.version$minor, sep = "."), cpu = cpu(),
         date = format(Sys.Date()))
ver <- function(ref) {
  p <- sub("::.*", "", ref)
  if (p == "base R") "" else format(utils::packageVersion(p))
}
if (length(only)) {
  cases <- Filter(function(k) k$fn %in% only, cases)
}
rows <- list()
for (k in cases) {
  for (n in k$n) {
    set.seed(20261112)
    e <- list2env(k$mk(n), parent = globalenv())
    row <- data.frame(fn = k$fn, case = k$case, ref = k$ref, n = n, of = k$of,
                      crimlab = NA_real_, other = NA_real_, ratio = NA_real_,
                      agree = "error", tol = k$tol,
                      scale = if (k$sc) "largest" else "each", refver = "")
    msg <- tryCatch({
      row$refver <- ver(k$ref)
      ra <- eval(k$a, e)
      rb <- eval(k$b, e)
      if (same(val(k$va, ra, e), val(k$vb, rb, e), k$tol, k$sc)) {
        row$agree <- "yes"
        tm <- bench::mark(exprs = list(k$a, k$b), env = e, check = FALSE,
                          min_iterations = 3, min_time = 0.5,
                          memory = FALSE, filter_gc = FALSE,
                          time_unit = "s")
        row$crimlab <- as.numeric(tm$median[1])
        row$other <- as.numeric(tm$median[2])
        row$ratio <- row$other / row$crimlab
      } else {
        row$agree <- "no"
      }
      ""
    }, error = function(x) conditionMessage(x))
    rows[[length(rows) + 1L]] <- row
    cat(sprintf("%-8s %-42s %8g %10.3g %10.3g %8.3g %s %s\n", k$fn, k$ref, n,
                row$crimlab, row$other, row$ratio, row$agree, msg))
  }
}
tab <- cbind(do.call(rbind, rows), t(env))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, file.path(out, "bench.csv"), row.names = FALSE)

g <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "g"))
md <- c(
  "<!-- written by dev/bench.R from bench.csv; do not edit -->",
  "",
  sprintf("R %s, %s, %s. Median seconds per call; ratio = other / crimlab (above 1: crimlab faster).",
          env[["r"]], env[["cpu"]], env[["date"]]),
  "",
  "| function | counterpart | n | of | crimlab | other | ratio | agree |",
  "|---|---|--:|---|--:|--:|--:|---|",
  sprintf("| %s | %s%s | %s | %s | %s | %s | %s | %s |",
          ifelse(nzchar(tab$case), paste0(tab$fn, " (", tab$case, ")"), tab$fn),
          tab$ref,
          ifelse(nzchar(tab$refver), paste0(" ", tab$refver), ""),
          format(tab$n, scientific = FALSE, big.mark = ","), tab$of,
          g(tab$crimlab), g(tab$other), g(tab$ratio), tab$agree)
)
writeLines(md, file.path(out, "bench.md"))

# the README lists, from the table: no figures
if (!length(only) && file.exists("README.md")) {
  pk <- tab$ref != "base R"
  fns <- unique(tab$fn)
  # faster means by at least a tenth: a smaller margin is within what two
  # runs of the same call differ by
  slow <- unique(tab$fn[!is.na(tab$ratio) & tab$ratio < 1.1])
  cp <- unique(tab$fn[pk])
  lst <- function(x) if (length(x)) paste0("`", x, "()`", collapse = ", ") else "none"
  bo <- setdiff(fns, cp)
  txt <- c(
    "<!-- bench start: written by dev/bench.R -->",
    paste0("Faster, by at least a tenth, than every package that computes ",
           "the same quantity, at every size timed: ", lst(setdiff(cp, slow)),
           "."),
    "",
    paste0("Not faster than such a package at some size timed: ",
           lst(intersect(cp, slow)), "."),
    "",
    paste0("No package found that computes the same quantity, so timed ",
           "against a base-R implementation of the formula: ", lst(bo),
           "; of these, not faster than base R at some size: ",
           lst(intersect(bo, slow)), "."),
    "<!-- bench end -->"
  )
  rd <- readLines("README.md")
  i <- grep("<!-- bench start", rd, fixed = TRUE)
  j <- grep("<!-- bench end -->", rd, fixed = TRUE)
  if (length(i) == 1L && length(j) == 1L && j > i) {
    writeLines(c(rd[seq_len(i - 1L)], txt, rd[-seq_len(j)]), "README.md")
  }
}
bad <- tab[tab$agree != "yes", c("fn", "ref", "n", "agree")]
if (nrow(bad)) {
  print(bad, row.names = FALSE)
  quit(status = 1L)
}
