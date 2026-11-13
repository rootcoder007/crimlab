# Speed of every exported function against every package found that
# computes the same quantity (dev/cpart.R records the search), on fixed
# seeded inputs. A pair is timed only after the two results agree to the
# stated relative tolerance; a pair that disagrees is reported as such and
# fails the run. Where no package computes the quantity, the counterpart is
# a plain base-R implementation of the same formula ("base R"). Calls are
# timed as a user makes them, argument checks included; graphics go to a
# null device. One process, no parallel workers.
#
#   Rscript dev/bench.R OUTDIR [fn ...]       one run, into OUTDIR/bench.csv
#   Rscript dev/bench.R --merge RUNDIR ...    the runs merged into
#                                             inst/bench/bench.csv, bench.md
#                                             and the README lists
#
# Each run records the machine, CPU, R version, CPU frequency governor
# (where the system reports one), the compiler flags R builds packages
# with, the 1-minute load average at its start and end, and the time; each
# row records the load before and after its timing, and a row whose load
# rose by more than 1.5 while it was timed is timed again (up to three
# times; the row keeps its last loads). dev/benchrun.sh takes the runs.
# The merge keeps, for each row, the run whose ratio is the median over the
# runs, with the smallest and largest ratio (lo, hi). Exit status 1 on a
# disagreement or an error, never on a slow result.
#
# Columns: fn, case (which arguments, where one function is timed two
# ways), ref (package::function), kind ("third" for another package, "own"
# for rmorie, the maintainer's own package, "base" for base R), n and of
# (what n counts), crimlab and other (median seconds per call), ratio
# (other / crimlab: above 1, crimlab is faster), agree ("yes", "no",
# "error", or "not run" with the reason in note), tol and scale (relative to
# each value or to the largest), refver (the counterpart's version), and the
# row's l0 and l1 (load before and after its timing), and the run's
# mach, cpu, rver, gov, load0, load1 and date (when it started), named as
# in the other speed tables (dev/bench-sp.R).
#
# Counterparts: rmorie (r-universe, not CRAN); Familias, forrel, paramlink
# and relSim for kinidx and patidx (Familias also gtmp); forensim and
# relSim for mixlr; comparison and idiolect (functions it bundles from the
# ROC package, unexported) for cllr, cllrmin and lrcal, comparison for
# eceplt; DasGuptR and DemoDecomp for clrsmp. forensim declares tkrplot,
# which needs Tk headers and, to load, an X display; where it cannot be
# loaded its pair is reported "not run".
a0 <- commandArgs(TRUE)
merge <- identical(a0[1], "--merge")
out <- if (length(a0)) a0[1] else "bench-run"
only <- a0[-1]
if (!merge) {
suppressMessages(library(crimlab))
grDevices::pdf(NULL)
}
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
  pick("shday", "shnight", "ratio"),
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
  pick("rtoff", "rtpair", "nulloff", "ratio"),
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
  # relSim's forms: a profile is the 2 x loci allele codes, the smaller
  # first at each locus (as relSim's own profiles are: it reads an unsorted
  # pair wrongly), the frequencies a list of one vector a locus
  rp <- function(g) as.integer(apply(g, 1, sort))
  list(fq = fq, g1 = g1, g2 = g2, f1 = fr(g1), mo = mo, fa = fa, ch = ch,
       fc = fr(ch), ff = fr(fa), r1 = rp(g1), r2 = rp(g2), rch = rp(ch),
       rfa = rp(fa), rf = list(freqs = lapply(seq_len(n), function(i) fq[i, ])))
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
# paramlink, forrel's predecessor: the same pedigrees, numbered
pl <- function(fq, h1, h2, gs) {
  ann <- lapply(seq_len(nrow(fq)), function(i) list(alleles = 1:6, afreq = fq[i, ]))
  set1 <- function(x) {
    ids <- x$orig.ids
    m <- matrix(0, length(ids), 2 * nrow(fq))
    for (id in names(gs)) if (id %in% ids) m[match(id, ids), ] <- c(t(gs[[id]]))
    paramlink::setMarkers(x, m, annotations = ann)
  }
  r <- paramlink::LR(list(set1(h1), lapply(h2, set1)), ref = 2)
  l <- r[["likelihoodsPerSystem"]]
  sum(log10(l[, 1] / l[, 2]))
}
plki <- function(fq, g1, g2) {
  pl(fq, paramlink::nuclearPed(2),
     list(paramlink::singleton(3), paramlink::singleton(4)),
     list("3" = g1, "4" = g2))
}
plpi <- function(fq, ch, fa, mo) {
  pl(fq, paramlink::nuclearPed(1),
     list(paramlink::relabel(paramlink::nuclearPed(1), c(4, 2, 3)),
          paramlink::singleton(1)),
     list("1" = fa, "2" = mo, "3" = ch))
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
add("kinidx", "paramlink::LR", c(20, 200), "loci", dna,
  kinidx(g1, g2, f1, "sibling"), plki(fq, g1, g2),
  function(r) r[["log10idx"]], identity)
add("patidx", "paramlink::LR", c(20, 200), "loci", dna,
  patidx(ch, fa, fc, mo), plpi(fq, ch, fa, mo),
  function(r) r[["log10idx"]], identity)
add("kinidx", "relSim::lrSib", c(20, 200), "loci", dna,
  kinidx(g1, g2, f1, "sibling"), relSim::lrSib(r1, r2, rf),
  function(r) r[["log10idx"]], log10)
# the child shares an allele with the alleged father at every locus
add("kinidx", "relSim::lrPC", c(20, 200), "loci", dna,
  kinidx(fa, ch, ff, "parent"), relSim::lrPC(rfa, rch, rf),
  function(r) r[["log10idx"]], log10, case = "parent")
add("patidx", "relSim::lrPC", c(20, 200), "loci", dna,
  patidx(ch, fa, fc), relSim::lrPC(rfa, rch, rf),
  function(r) r[["log10idx"]], log10, case = "no mother")

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

# victim (alleles 1, 2) and suspect (3, 4) against the victim and one
# unknown, one locus
add("mixlr", "relSim::lrMix", 4, "alleles",
  function(n) {
    f <- stats::runif(n, 0.02, 0.08)
    list(p = stats::setNames(f, seq_len(n)),
         pr = list(list(prof1 = 1:2, prof2 = 3:4)),
         fq = list(freqs = list(c(f, 1 - sum(f))), loci = "L1"))
  },
  mixlr(p, character(0), 0, c("3", "4"), 1), relSim::lrMix(pr, fq),
  function(r) r[["lr"]], c, case = "two people")

lrs <- function(n) {
  s <- exp(stats::rnorm(n / 2, 1, 1.5))
  d <- exp(stats::rnorm(n / 2, -1, 1.5))
  # idiolect's form: natural-log scores with the truth of each comparison
  list(s = s, d = d, x = data.frame(score = log(c(s, d)),
                                    target = rep(c(TRUE, FALSE), each = n / 2)))
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
add("cllr", "idiolect:::Cllr", c(1e3, 1e5, 1e6), "ratios", lrs,
  cllr(s, d), idiolect:::Cllr(x), identity, identity)
add("cllrmin", "idiolect:::Cllr", c(1e3, 1e5, 1e6), "ratios", lrs,
  cllrmin(s, d), idiolect:::Cllr(x, TRUE), pick("min"), identity)
add("lrcal", "idiolect:::opt.llr", c(1e3, 1e5, 1e6), "ratios", lrs,
  lrcal(s, d), idiolect:::opt.llr(x, FALSE),
  pick("same", "diff"), function(r) exp(r[["opt.llr"]]))
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

# DasGuptR: the rate and the share of each part as two factors of
# sum(rate * share); DemoDecomp::horiuchi on the same function, exact for
# it, whose cost grows with the square of the parts
dgx <- function(n) {
  ra <- stats::runif(n, 10, 1e3)
  rb <- stats::runif(n, 10, 1e3)
  ca <- ra * stats::runif(n)
  cb <- rb * stats::runif(n)
  list(ra = ra, rb = rb, ca = ca, cb = cb,
       x = data.frame(pop = rep(c("X1", "X2"), each = n), k = rep(seq_len(n), 2),
                      clr = c(ca / ra, cb / rb), shr = c(ra / sum(ra), rb / sum(rb))),
       p1 = c(ca / ra, ra / sum(ra)), p2 = c(cb / rb, rb / sum(rb)),
       f = function(p) {
         k <- length(p) / 2
         sum(p[seq_len(k)] * p[k + seq_len(k)])
       })
}
add("clrsmp", "DasGuptR::dgnpop", c(5, 1e3, 1e5), "offence types", dgx,
  clrsmp(ca, ra, cb, rb),
  suppressMessages(DasGuptR::dgnpop(x, pop = "pop", factors = c("clr", "shr"),
                                    id_vars = "k", ratefunction = "sum(clr*shr)")),
  function(r) r$total[c("rateeff", "mixeff")],
  function(r) {
    v <- function(f, p) r$rate[r$factor == f & r$pop == p]
    c(v("clr", "X2") - v("clr", "X1"), v("shr", "X2") - v("shr", "X1"))
  }, sc = TRUE)
add("clrsmp", "DemoDecomp::horiuchi", c(5, 1e3), "offence types", dgx,
  clrsmp(ca, ra, cb, rb), DemoDecomp::horiuchi(f, p1, p2, N = 20),
  function(r) c(r$parts$rateeff, r$total[["mixeff"]]),
  function(r) {
    k <- length(r) / 2
    c(r[seq_len(k)], sum(r[k + seq_len(k)]))
  }, sc = TRUE)

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
  sntbnd(y, z, "custody"), rmorie::morie_sentence_effect_bounds(y, z),
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
sh <- function(cmd, args = character()) {
  tryCatch(suppressWarnings(system2(cmd, args, stdout = TRUE, stderr = FALSE)),
           error = function(e) character(0))
}
cpu <- function() {
  m <- sub("^Model name: *", "", grep("^Model name:", sh("lscpu"), value = TRUE))
  if (!length(m)) m <- sh("sysctl", c("-n", "machdep.cpu.brand_string"))
  if (length(m)) paste(m[1], Sys.info()[["machine"]]) else Sys.info()[["machine"]]
}
rd1 <- function(f) {
  x <- tryCatch(suppressWarnings(readLines(f, n = 1L)), error = function(e) "")
  if (length(x)) x else ""
}
# the 1-minute load average, as uptime reports it
load1 <- function() {
  x <- rd1("/proc/loadavg")
  if (!nzchar(x)) x <- gsub("[{}]", "", paste(sh("sysctl", c("-n", "vm.loadavg")), collapse = " "))
  as.numeric(strsplit(trimws(x), " +")[[1]][1])
}
kind <- function(ref) {
  ifelse(ref == "base R", "base", ifelse(startsWith(ref, "rmorie::"), "own", "third"))
}
pkg <- function(ref) sub(":::?.*", "", ref)

if (!merge) {
  env <- c(mach = Sys.info()[["nodename"]], cpu = cpu(),
           rver = R.version.string,
           gov = rd1("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor"),
           load0 = format(load1()),
           date = format(Sys.time(), "%Y-%m-%d %H:%M %Z"))
  flags <- sh(file.path(R.home("bin"), "R"), c("CMD", "config", "CXX17FLAGS"))
  # the builds timed
  ps <- unique(c("crimlab", setdiff(pkg(vapply(cases, `[[`, "", "ref")), "base R")))
  for (p in ps) {
    cat(sprintf("%-12s %s\n", p, tryCatch(paste(find.package(p),
      utils::packageVersion(p)), error = function(e) "not installed")))
  }
  if (length(only)) cases <- Filter(function(k) k$fn %in% only, cases)
  rows <- list()
  for (k in cases) {
    for (n in k$n) {
      set.seed(20261112)
      row <- data.frame(fn = k$fn, case = k$case, ref = k$ref, kind = kind(k$ref),
                        n = n, of = k$of, crimlab = NA_real_, other = NA_real_,
                        ratio = NA_real_, agree = "error", tol = k$tol,
                        scale = if (k$sc) "largest" else "each", refver = "",
                        note = "", l0 = NA_real_, l1 = NA_real_)
      p <- pkg(k$ref)
      if (p != "base R" && !have(p)) {
        row$agree <- "not run"
        row$note <- paste(p, "is not installed on this machine")
        msg <- ""
      } else {
        msg <- tryCatch({
          e <- list2env(k$mk(n), parent = globalenv())
          if (p != "base R") row$refver <- format(utils::packageVersion(p))
          ra <- eval(k$a, e)
          rb <- eval(k$b, e)
          if (same(val(k$va, ra, e), val(k$vb, rb, e), k$tol, k$sc)) {
            row$agree <- "yes"
            for (try in 1:3) {
              row$l0 <- load1()
              tm <- bench::mark(exprs = list(k$a, k$b), env = e, check = FALSE,
                                min_iterations = 3, min_time = 0.5,
                                memory = FALSE, filter_gc = FALSE,
                                time_unit = "s")
              row$l1 <- load1()
              if (!(row$l1 - row$l0 > 1.5)) break
            }
            row$crimlab <- as.numeric(tm$median[1])
            row$other <- as.numeric(tm$median[2])
            row$ratio <- row$other / row$crimlab
          } else {
            row$agree <- "no"
          }
          ""
        }, error = function(x) conditionMessage(x))
        if (nzchar(msg)) row$note <- msg
      }
      rows[[length(rows) + 1L]] <- row
      cat(sprintf("%-8s %-36s %8g %10.3g %10.3g %8.3g %s %s\n", k$fn, k$ref, n,
                  row$crimlab, row$other, row$ratio, row$agree, row$note))
    }
  }
  env[["load1"]] <- format(load1())
  tab <- cbind(do.call(rbind, rows), t(env))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(tab, file.path(out, "bench.csv"), row.names = FALSE)
  writeLines(flags, file.path(out, "flags.txt"))
  bad <- tab[tab$agree %in% c("no", "error"), c("fn", "ref", "n", "agree", "note")]
  if (nrow(bad)) {
    print(bad, row.names = FALSE)
    quit(status = 1L)
  }
  quit(status = 0L)
}

## merge ------------------------------------------------------------------

runs <- a0[-1]
rs <- lapply(runs, function(d) utils::read.csv(file.path(d, "bench.csv"),
                                                stringsAsFactors = FALSE))
key <- function(t) paste(t$fn, t$case, t$ref, t$n, sep = "\r")
k0 <- key(rs[[1]])
for (t in rs[-1]) {
  if (!identical(key(t), k0)) stop("the runs do not hold the same rows")
}
rat <- sapply(rs, function(t) t$ratio)
if (is.null(dim(rat))) rat <- matrix(rat, ncol = length(rs))
tab <- rs[[1]]
for (i in seq_len(nrow(tab))) {
  r <- rat[i, ]
  if (all(is.na(r))) next
  # the run holding the median ratio (the lower middle one for an even count)
  o <- order(r)
  j <- o[ceiling(sum(!is.na(r)) / 2)]
  tab[i, ] <- rs[[j]][i, ]
}
tab$lo <- suppressWarnings(apply(rat, 1, min, na.rm = TRUE))
tab$hi <- suppressWarnings(apply(rat, 1, max, na.rm = TRUE))
tab$lo[!is.finite(tab$lo)] <- NA
tab$hi[!is.finite(tab$hi)] <- NA
tab$runs <- rowSums(!is.na(rat))
dir.create("inst/bench", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, "inst/bench/bench.csv", row.names = FALSE)

g <- function(x) ifelse(is.na(x), "", formatC(x, digits = 3, format = "g"))
st <- vapply(rs, function(t) {
  sprintf("%s, load %s to %s", t$date[1], t$load0[1], t$load1[1])
}, "")
flags <- readLines(file.path(runs[1], "flags.txt"))
md <- c(
  "<!-- written by dev/bench.R from bench.csv; do not edit -->",
  "",
  sprintf("Machine %s: %s, %s%s; crimlab built with R's default flags (%s).",
          tab$mach[1], tab$cpu[1], tab$rver[1],
          if (!is.na(tab$gov[1]) && nzchar(tab$gov[1])) paste0(", CPU frequency governor ", tab$gov[1]) else "",
          paste(flags, collapse = " ")),
  "",
  sprintf("%d runs (start time, 1-minute load average at the start and at the end): %s. A row whose load rose by more than 1.5 while it was timed was timed again; the loads of each row are in bench.csv (l0, l1). Each row is from the run whose ratio is the median of the %d; lo and hi are the smallest and largest ratio over the runs.%s",
          length(rs), paste(st, collapse = "; "), length(rs),
          if (nzchar(Sys.getenv("BENCH_NOTE"))) paste0(" ", Sys.getenv("BENCH_NOTE")) else ""),
  "",
  "Median seconds per call; ratio = other / crimlab (above 1: crimlab faster). Kind: third, another package; own, rmorie (the maintainer's own package); base, a base-R implementation of the formula where no package computes the quantity.",
  "",
  "| function | counterpart | kind | n | of | crimlab | other | ratio | lo | hi | agree |",
  "|---|---|---|--:|---|--:|--:|--:|--:|--:|---|",
  sprintf("| %s | %s%s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
          ifelse(nzchar(tab$case), paste0(tab$fn, " (", tab$case, ")"), tab$fn),
          tab$ref, ifelse(nzchar(tab$refver), paste0(" ", tab$refver), ""),
          tab$kind, format(tab$n, scientific = FALSE, big.mark = ","), tab$of,
          g(tab$crimlab), g(tab$other), g(tab$ratio), g(tab$lo), g(tab$hi),
          ifelse(nzchar(tab$note), paste0(tab$agree, ": ", tab$note), tab$agree))
)
writeLines(md, "inst/bench/bench.md")

# the README lists, from the table: no figures. Faster means by at least a
# tenth (ratio 1.1) in the median run; a smaller margin is within what two
# runs of the same call differ by.
if (file.exists("README.md")) {
  lst <- function(x) if (length(x)) paste0("`", x, "()`", collapse = ", ") else "none"
  fns <- unique(tab$fn)
  ok <- tab$agree == "yes"
  slow <- function(kd) unique(tab$fn[ok & tab$kind == kd & tab$ratio < 1.1])
  has <- function(kd) unique(tab$fn[ok & tab$kind == kd])
  th <- has("third")
  ow <- setdiff(has("own"), th)
  bo <- setdiff(fns, c(th, ow))
  txt <- c(
    "<!-- bench start: written by dev/bench.R -->",
    paste0("Faster, by at least a tenth, than every other package found that ",
           "computes the same quantity, at every size timed: ",
           lst(setdiff(th, slow("third"))), "."),
    "",
    paste0("Not faster than another package at some size timed: ",
           lst(intersect(th, slow("third"))), "."),
    "",
    paste0("No other package found that computes the same quantity, so no ",
           "claim of speed against the field: ", lst(c(ow, bo)), ". Of these, ",
           "compared with rmorie, the maintainer's own package: faster at ",
           "every size timed: ", lst(setdiff(ow, slow("own"))),
           "; not faster at some size: ", lst(intersect(ow, slow("own"))),
           ". Compared with a base-R implementation of the formula: faster ",
           "at every size timed: ", lst(setdiff(bo, slow("base"))),
           "; not faster at some size: ", lst(intersect(bo, slow("base"))), "."),
    "<!-- bench end -->"
  )
  rd <- readLines("README.md")
  i <- grep("<!-- bench start", rd, fixed = TRUE)
  j <- grep("<!-- bench end -->", rd, fixed = TRUE)
  if (length(i) == 1L && length(j) == 1L && j > i) {
    writeLines(c(rd[seq_len(i - 1L)], txt, rd[-seq_len(j)]), "README.md")
  }
}
