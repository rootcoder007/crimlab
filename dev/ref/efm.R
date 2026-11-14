# Reference values from EuroForMix (github.com/oyvble/euroformix, the
# package euroformix) for tests/testthat/test-gmixlr.R, written to
# tests/testthat/ref/efm.txt: on small cases, the log likelihood of each
# account at given parameters (the compiled likelihood of prepareCobj(),
# with no restriction of the genotypes) and its maximum by calcMLE()
# (resttol = 0, normalize = FALSE, no forward stutter). Each case holds
# its inputs as gmixlr() takes them. The commit of EuroForMix is taken
# from EFM_COMMIT. From the package root, with euroformix in the library:
#   EFM_COMMIT=<sha> Rscript dev/ref/efm.R
# the packages it needs (dev/ref/all.R skips it without them)
need <- "euroformix"
source(file.path("dev", "ref", "put.R"))
suppressPackageStartupMessages(library(euroformix))
kit <- "ESX17"
ki <- getKit(kit)

# a case: peaks and references as data frames, frequencies by marker
mk <- function(peaks, freq, ref, hp, hd, nc, pc = 0.05, lambda = 0.01, theta = 0,
               stut = FALSE, deg = FALSE, par, fit = TRUE) {
  at <- 50
  loc <- unique(freq$locus)
  samples <- list(s1 = lapply(stats::setNames(loc, toupper(loc)), function(l) {
    p <- peaks[peaks$locus == l, ]
    list(adata = as.character(p$allele), hdata = p$height)
  }))
  pf <- lapply(stats::setNames(loc, toupper(loc)), function(l) {
    f <- freq[freq$locus == l, ]
    stats::setNames(f$freq, as.character(f$allele))
  })
  ids <- unique(ref$id)
  rd <- lapply(stats::setNames(ids, ids), function(i) {
    r <- ref[ref$id == i, ]
    lapply(stats::setNames(seq_len(nrow(r)), toupper(r$locus)), function(k) {
      as.character(c(r$a1[k], r$a2[k]))
    })
  })
  if (deg) {
    kt <- ki[toupper(ki$Marker) %in% toupper(loc), ]
    freq$size <- kt$Size[match(paste(toupper(freq$locus), freq$allele),
                               paste(toupper(kt$Marker), kt$Allele))]
  }
  one <- function(who) {
    co <- match(ids, who, nomatch = 0)
    kn <- which(co == 0)
    if (!length(kn)) kn <- NULL
    cc <- euroformix:::prepareC(nc, samples, pf, rd, co, kn, if (deg) kit else NULL,
                                stut, FALSE, at, pc, lambda, theta, NULL, NULL, NULL,
                                FALSE, FALSE)
    cc$maxThreads <- 1L
    cc$nU <- nc - length(who)
    obj <- prepareCobj(cc)
    th <- c(par[paste0("m", seq_len(nc))], par["mu"], par["omega"],
            if (deg) par["beta"] else 1, if (stut) par["xi"] else 0, 0)
    obj$calcGenoWeightsMax(as.numeric(th))
    obj$restrictGenos(0)
    ll <- obj$loglik(as.numeric(th))
    obj$close()
    mle <- NA_real_
    if (fit) {
      m <- calcMLE(nc, samples, pf, rd, co, kn, if (deg) kit else NULL, DEG = deg, BWS = stut,
                   FWS = FALSE, AT = at, pC = pc, lambda = lambda, fst = theta, normalize = FALSE,
                   resttol = 0, seed = 1, nDone = 3)
      mle <- m$fit$loglik
    }
    c(ll = ll, mle = mle)
  }
  h1 <- one(hp)
  h2 <- one(hd)
  list(peaks = peaks, freq = freq, ref = ref, hp = hp, hd = hd, nc = nc, pc = pc,
       lambda = lambda, theta = theta, stut = stut, deg = deg, par = par,
       ll = unname(c(h1["ll"], h2["ll"])), mle = unname(c(h1["mle"], h2["mle"])))
}

# frequencies of two markers of the kit, summing to 1 at each, the largest
# allele of the kit at each marker listed (EuroForMix sizes the alleles
# not seen as that allele)
f2 <- data.frame(locus = rep(c("D3S1358", "VWA"), c(6, 6)),
                 allele = c(13, 14, 15, 16, 17, 21, 15, 16, 17, 18, 19, 25),
                 freq = c(0.1, 0.2, 0.25, 0.25, 0.15, 0.05, 0.1, 0.25, 0.3, 0.2, 0.1, 0.05))
r2 <- data.frame(id = rep(c("V", "S"), each = 2), locus = rep(c("D3S1358", "VWA"), 2),
                 a1 = c(14, 16, 15, 17), a2 = c(16, 17, 17, 18))
pk <- data.frame(locus = rep(c("D3S1358", "VWA"), c(4, 3)), allele = c(14, 15, 16, 17, 16, 17, 18),
                 height = c(1210, 380, 1150, 420, 1400, 1960, 610))
cases <- list(
  # one contributor, both peaks, no drop-in
  one = mk(pk[c(1, 3), ], f2, r2[r2$id == "V", ], "V", character(0), 1, pc = 0,
           par = c(mu = 1100, omega = 0.15, m1 = 1)),
  # one contributor with a drop-in peak and a dropped-out allele
  din = mk(data.frame(locus = c("D3S1358", "D3S1358", "VWA"), allele = c(15, 13, 17),
                      height = c(900, 120, 1000)),
           f2, r2[r2$id == "S", ], "S", character(0), 1,
           par = c(mu = 1000, omega = 0.25, m1 = 1)),
  # two contributors, the victim known under both accounts
  two = mk(pk, f2, r2, c("V", "S"), "V", 2, par = c(mu = 1800, omega = 0.12, m1 = 0.73, m2 = 0.27)),
  # the minor contributor's 18 at vWA dropped out
  drop = mk(pk[-7, ], f2, r2, c("V", "S"), "V", 2, par = c(mu = 1700, omega = 0.2, m1 = 0.8, m2 = 0.2)),
  # backward stutter peaks one repeat below the major alleles
  stut = mk(rbind(pk, data.frame(locus = c("D3S1358", "VWA"), allele = c(13, 15), height = c(95, 110))),
            f2, r2, c("V", "S"), "V", 2, stut = TRUE,
            par = c(mu = 1800, omega = 0.12, xi = 0.07, m1 = 0.73, m2 = 0.27)),
  # degradation by fragment size (the kit's sizes), with stutter
  deg = mk(rbind(pk, data.frame(locus = "D3S1358", allele = 13, height = 95)), f2, r2, c("V", "S"),
           "V", 2, stut = TRUE, deg = TRUE,
           par = c(mu = 1900, omega = 0.13, beta = 0.8, xi = 0.06, m1 = 0.7, m2 = 0.3)),
  # coancestry
  fst = mk(pk, f2, r2, c("V", "S"), "V", 2, theta = 0.03,
           par = c(mu = 1800, omega = 0.12, m1 = 0.73, m2 = 0.27)),
  # two unknowns under the second account
  unk2 = mk(pk, f2, r2, c("V", "S"), character(0), 2, fit = FALSE,
            par = c(mu = 1800, omega = 0.15, m1 = 0.7, m2 = 0.3))
)
cases$commit <- Sys.getenv("EFM_COMMIT", "unknown")
put(cases, "efm", need)
