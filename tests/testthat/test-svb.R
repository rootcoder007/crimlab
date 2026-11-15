# Bayesian scaling (src/svb.cpp): chains drawn from the same seed equal the
# reference chains of the rmorie samplers they port. A chain started from an
# eigenvector or singular vector equals the reference one when that vector's
# sign (arbitrary in both) agrees, as it does on these inputs. Each draw of a
# chain depends on the one before, so rounding differences carry along the
# chain; the cjrirt() chain is compared at 1e-9, the rest at 1e-10.

test_that("bamscl() equals the reference chain", {
  set.seed(2)
  tr <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(vapply(1:12, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * tr + rnorm(5, 0, 0.3)
  }, numeric(5)))
  r <- bamscl(Z, 10, 10)
  # rmorie 1.5.0: R$morie_spatial_voting_bayesian_am(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-1.3717372738145706, -0.456070460967687, -0.058626758592582533,
        0.70485113534015764, 1.1815833580346824)
  expect_equal(unname(c(r$zeta)), c(ref), tolerance = 1e-10)
  ref <- c(0.032351879017836101, 0.047714116677280309, 0.089602057876970412,
        0.050902133238955108, 0.039392345664905203)
  expect_equal(unname(c(r$sd)), c(ref), tolerance = 1e-10)
  ref <- c(4.1503175442802913, 4.0884914045526113, 3.8858491494199092,
        3.9831703984429234, 4.7093761489547639, 3.80856739792373,
        4.0493814674599955, 3.8900888664869342, 3.9723453732391363,
        3.9460213273298828, 3.9299064624508366, 4.0542989443536053)
  expect_equal(unname(c(r$a)), c(ref), tolerance = 1e-10)
  ref <- c(1.607959412484447, 14.134467533074556, 7.675513142305288,
        13.341988264246616, 4.3969824900691723)
  expect_equal(unname(c(r$tau)), c(ref), tolerance = 1e-10)
})


test_that("bmdsf() equals the reference chain", {
  set.seed(7)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  r <- bmdsf(D, 2, 5, 5)
  # rmorie 1.5.0: R$morie_spatial_voting_bayesian_mds(), the same input,
  # its values printed with 17 significant digits
  ref <- c(2.071720377634779, -1.1849358962043985, -0.75500315045888933,
        -0.97979711246753065, -0.88474299571158677, -1.1988404886714361,
        -2.787032299189101, -1.8446581724321072, -1.3261647971872881,
        0.43402461715217711, -1.3635280183332703, 0.93698326597392045)
  expect_equal(unname(c(r$pos)), c(ref), tolerance = 1e-10)
  ref <- c(0, 3.4246674128304067, 3.2348968655973378, 4.4568518584753418,
        3.3108632125010047, 4.9659835663537137, 3.4246674128304067, 0,
        0.74310444271111264, 2.3901388327557149, 0.60086582614917139,
        2.8030694993230374, 3.2348968655973378, 0.74310444271111264, 0,
        1.853213597861965, 0.31600270672791403, 2.3637266028315329,
        4.4568518584753418, 2.3901388327557149, 1.853213597861965, 0,
        1.9019387118786448, 1.0122807145171038, 3.3108632125010047,
        0.60086582614917139, 0.31600270672791403, 1.9019387118786448, 0,
        2.3604658762297985, 4.9659835663537137, 2.8030694993230374,
        2.3637266028315329, 1.0122807145171038, 2.3604658762297985, 0)
  expect_equal(unname(c(r$dist)), c(ref), tolerance = 1e-10)
  ref <- 0.34050019837511436
  expect_equal(unname(c(r$sigma)), c(ref), tolerance = 1e-10)
  ref <- 8.6975785872760909
  expect_equal(unname(c(r$tau)), c(ref), tolerance = 1e-10)
})


test_that("bunfld() equals the reference chain", {
  set.seed(2)
  S <- matrix(runif(6, -0.5, 0.5), 3)
  X <- matrix(runif(16, -0.5, 0.5), 8)
  D <- sqrt(outer(X[, 1], S[, 1], "-")^2 + outer(X[, 2], S[, 2], "-")^2)
  r <- bunfld(D, 2, 5, 5)
  # rmorie 1.5.0: R$morie_spatial_voting_bayesian_unfolding(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-1.3162085014919538, -1.5519129052647065, -2.2937746672884272,
        0.081275419307444557, 0.84998229951063475, 0.1673992766492842)
  expect_equal(unname(c(r$stim)), c(ref), tolerance = 1e-10)
  ref <- c(-1.2897221482803696, -1.8558840189220021, -2.4533838649828597,
        -1.2487210974783627, -1.1772068122135217, -1.5664459876016772,
        -2.0379755279462435, -1.2560153907751159, 0.43122199488585533,
        0.078145030839969043, 0.40933213905955107, 0.63578014798778948,
        -0.023821248500059754, -0.028794162604498651, 0.59920502734555225,
        0.074758023145452018)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
  ref <- c(0.57177534293879162, 1.1065447058909419, 1.2595779785980112,
        1.7861402115486378, 0.81709807930254319, 0.45708422835410351,
        1.3021951479208762, 0.72142482795081031, 0.86870212165989913,
        1.0190904691254663, 1.1300696211167722, 2.0642912737847139,
        1.0781095950190818, 1.1050059820552571, 0.9809830911698374,
        0.98158289340147908, 1.0848968165256441, 0.79821031417814226,
        0.45784594002615453, 1.8954111080184943, 1.2370140162318335,
        0.87479376731842784, 1.2492806453789069, 1.068089504055926)
  expect_equal(unname(c(r$dist)), c(ref), tolerance = 1e-10)
  ref <- 1.1736544751535549
  expect_equal(unname(c(r$tau)), c(ref), tolerance = 1e-10)
})


test_that("cjrirt() equals the reference Philox chain", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  r <- cjrirt(v, 1, 10, 5, seed = 3)
  # rmorie 1.5.0: R$morie_spatial_voting_cjr_irt(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-0.16313307857480586, 0.19561755676967202, -0.33365693856170353,
        -0.44426271212974466, -0.4229324910604102, -0.45435092044093028,
        0.15548270587770932, -0.04068017699158645)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-09)
  ref <- c(0.36738859490441911, 1.4581424381422794, -0.1092697192402767,
        1.3677208542676191, -1.4386416659036469, 0.24582709140715697,
        0.86020829077960714, -0.058944835169904966, 0.021290392387438802,
        1.2224953784580941)
  expect_equal(unname(c(r$alpha)), c(ref), tolerance = 1e-09)
  ref <- c(-5.4264888214351306, -5.6950777111259479, -3.8135224306888866,
        -2.8579455735486001, 0.82764463545031541, -0.077839227510720463,
        -1.7087522788800809, -1.8778377760127947, -6.62020339011627,
        -0.45210032586618953)
  expect_equal(unname(c(r$beta)), c(ref), tolerance = 1e-09)
})


test_that("anomin() equals the reference chain", {
  set.seed(12)
  x <- seq(-1, 1, length.out = 10)
  V <- outer(x, runif(16, -0.7, 0.7), ">") * 1
  fl <- sample(length(V), 16)
  V[fl] <- 1 - V[fl]
  r <- anomin(V, 1, 4, 2, minvotes = 8, polarity = 10)
  # rmorie 1.5.0: R$morie_spatial_voting_alpha_nominate(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-2.0020075241105615, -1.7311865460372666, -0.45653372000626524,
        -0.14219080207266854, 0.5772877109324337, 0.76982962254188547,
        0.39185698359319404, 0.51942951935304482, 1.5443963069476951,
        0.52911844885850867)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
  ref <- 0.7600354995875469
  expect_equal(unname(c(r$alpha)), c(ref), tolerance = 1e-10)
  ref <- 7.3053960698220504
  expect_equal(unname(c(r$beta)), c(ref), tolerance = 1e-10)
  ref <- c(0.19169672886556366, 0.25574057426113334, 1.9202412411771252,
        -1.8292990235459887, 0.77143676902865754, -0.68100816620362226,
        0.42830139770767928, 0.45855266311913234, 0.42258435966217073,
        0.050343150961000414, -0.20420926440061815, 0.40580331239969425,
        -1.0016848277814461, 0.56537876130161024, 0.42606010923642512,
        0.18930412728439727)
  expect_equal(unname(c(r$yea)), c(ref), tolerance = 1e-10)
})


test_that("ordirt() equals the reference chain", {
  set.seed(2)
  th <- rnorm(15)
  Y <- matrix(cut(outer(th, runif(4, 0.8, 1.5)) + rnorm(60),
    c(-Inf, -0.5, 0.5, Inf), labels = FALSE), 15)
  r <- ordirt(Y, 1, 5, 5)
  # rmorie 1.5.0: R$morie_spatial_voting_ordinal_irt(), the same input,
  # its values printed with 17 significant digits
  ref <- c(-0.20711094761972862, -0.071862771240624859, 0.99536186299794416,
        -0.65969234239864438, -0.26670311303461502, -0.24445693995480325,
        0.11666131645892101, -0.99137049573990021, 0.35952957937841434,
        -0.45154642726156735, 0.088315975326005047, 0.60517626720490958,
        0.14934138239786954, -0.1302671115472005, 1.1869119847464138)
  expect_equal(unname(c(r$ideal)), c(ref), tolerance = 1e-10)
  ref <- c(0.66990746775333676, 1.4494736297153135, 2.4564108956516808,
        1.5945928700837184)
  expect_equal(unname(c(r$disc)), c(ref), tolerance = 1e-10)
  ref <- c(0.63921568326131417, -0.017606217791777796, 0.81543772770073786,
        0.21937960161625084)
  expect_equal(unname(c(r$icept)), c(ref), tolerance = 1e-10)
  # one vector of cutpoints per item, as in rmorie
  ref <- c(0, 0.51320773329748204, 0, 0.48738314450348674, 0,
    0.46753967202491731, 0, 0.50924824205705455)
  expect_equal(unlist(r$cuts), ref, tolerance = 1e-10)
})


test_that("dynirt() equals the reference chain", {
  # The chain starts from the first principal component of the votes,
  # whose sign is arbitrary: crimlab makes its largest loading positive,
  # rmorie takes the sign LAPACK gives. On this input the two agree. (On
  # the input of set.seed(14) they do not, and the chains then differ from
  # the first sweep; rmorie's sampler started from the reflected component
  # gives crimlab's chain there to 1e-15.)
  set.seed(15)
  per <- rep(1:3, each = 4)
  V <- matrix(rbinom(72, 1, 0.5), 6)
  V[2, 5] <- NA
  rot <- prcomp(replace(V, is.na(V), 0.5))$rotation[, 1]
  skip_if(rot[which.max(abs(rot))] < 0, "LAPACK gives the other sign")
  r <- dynirt(V, per, 8, 4, c0 = 1, d0 = 1)
  # rmorie 1.5.0: R$morie_spatial_voting_dynamic_irt(), the same input,
  # its values printed with 17 significant digits
  ref <- c(0.38310770340618888, 0.70107336274899745, 0.25878558556496961,
    -0.17035657826539177, 0.96757437081786302, 0.5323069400141307,
    -0.60321299446323862, 0.50799940604980787, -0.021349635579193369,
    0.25868360230260079, 0.80365989519072012, -0.23484570066134691,
    -1.4256119992224869, 0.41912911283156534, -0.15618560158209968,
    -0.080210634495035252, -0.428578663652187, -0.68728676370188402)
  expect_equal(unname(c(r$theta)), ref, tolerance = 1e-10)
  ref <- c(-1.0121492607245173, -4.0016337133816524, -2.1086315746335722,
    -1.189019551023718, -0.20622898787958033, 0.15025395298632915,
    0.4567927241254508, 0.51275527171353441, -1.1836630579831739,
    -1.7549490107770502, -2.2381011822944874, -0.18867468907781204)
  expect_equal(unname(c(r$alpha)), ref, tolerance = 1e-10)
  ref <- c(-2.3025198281026422, -7.1225112599092677, -2.4839223505686263,
    -1.6585430384257482, -0.098049599593633643, 0.30559601063768504,
    2.3943971144116367, -1.9976014264182593, -2.8257189023912566,
    2.7574343152324454, 2.3991517320802913, -2.2996838460148532)
  expect_equal(unname(c(r$beta)), ref, tolerance = 1e-10)
  ref <- c(1.7446268661944062, 3.0195581205270954, 0.46021856058776561,
    1.0411010003484291, 1.9369779672829677, 1.5345057441925645)
  expect_equal(unname(c(r$tau2)), ref, tolerance = 1e-10)
})


# Behaviour of the Bayesian scaling samplers, carried over from the tests
# of the methods they port.

test_that("cjrirt() gives the Philox chain of the reference sampler", {
  V <- em_votes()
  f <- cjrirt(V, nsamp = 200, burn = 50, seed = 5)
  got <- c(f$ideal[1:3, 1], f$alpha[1:2], f$beta[1, 1], f$sd[1, 1])
  # the same chain drawn by an independent implementation of the sampler
  ref <- c(-0.577500384336, -0.678132861519, -0.176501975313, -1.12013879797,
           0.540891603986, -9.235856442311, 0.297621796656)
  expect_lt(max(abs(got - ref)), 1e-9)
  expect_identical(dim(f$xchain), c(200L, 30L, 1L))
  g <- cjrirt(V, nsamp = 20, burn = 5, seed = 5, chains = FALSE)
  expect_null(g$xchain)
  expect_error(cjrirt(V, seed = 1.5), "whole number")
})

test_that("bamscl() recovers standardised stimulus positions", {
  set.seed(1)
  zt <- c(-1, -0.3, 0.4, 1.2)
  Z <- rnorm(60, 0, 0.5) + outer(runif(60, 0.5, 1.5), zt) +
    matrix(rnorm(240, 0, 0.3), 60, 4)
  f <- bamscl(Z, nsamp = 100, burn = 100)
  expect_gt(cor(f$zeta, zt), 0.99)
  expect_lt(f$zeta[1], 0)
  expect_identical(dim(f$draws), c(100L, 4L))
  expect_true(all(f$interval[1, ] <= f$interval[2, ]))
  expect_error(bamscl(matrix(rnorm(10), 10)), "at least two stimuli")
  expect_error(bamscl(cbind(NA, Z)), "no placements")
})

test_that("bmdsf() recovers the distances and aligns its draws", {
  set.seed(2)
  D <- as.matrix(dist(matrix(rnorm(16), 8)))
  f <- bmdsf(D, nsamp = 60, burn = 60)
  expect_gt(cor(f$dist[upper.tri(D)], D[upper.tri(D)]), 0.95)
  expect_identical(dim(f$draws), c(60L, 8L, 2L))
  expect_equal(f$pos, apply(f$draws, c(2, 3), mean),
               tolerance = testthat_tolerance())
  expect_error(bmdsf(matrix(0, 3, 3)), "positive distances")
  expect_error(bmdsf(matrix(1, 2, 3)), "square")
})

test_that("bunfld() recovers the respondent-stimulus geometry", {
  set.seed(3)
  S <- matrix(runif(8, -0.5, 0.5), 4)
  R <- matrix(runif(60, -0.5, 0.5), 30)
  D <- sqrt(outer(R[, 1], S[, 1], "-")^2 + outer(R[, 2], S[, 2], "-")^2)
  f <- bunfld(D, nsamp = 60, burn = 60)
  expect_gt(cor(as.numeric(f$dist), as.numeric(D)), 0.9)
  expect_identical(dim(f$ideal), c(30L, 2L))
  expect_error(bunfld(matrix(0, 3, 3)), "positive dissimilarities")
})

test_that("anomin() keeps the order of simulated legislators", {
  set.seed(4)
  x <- seq(-1, 1, length.out = 20)
  V <- outer(x, runif(40, -0.7, 0.7), ">") * 1
  flip <- sample(length(V), 40)
  V[flip] <- 1 - V[flip]
  f <- anomin(V, nsamp = 30, burn = 20, minvotes = 10, polarity = 20)
  expect_gt(cor(f$ideal[, 1], x[f$legis]), 0.95)
  expect_gte(f$alpha, 0)
  expect_lte(f$alpha, 1)
  g <- anomin(V, nsamp = 5, burn = 2, minvotes = 10, polarity = 20,
              fixalpha = TRUE)
  expect_identical(unique(g$draws$alpha), 1)
  expect_error(anomin(V * 2), "1 \\(yea\\), 0 \\(nay\\)")
  expect_error(anomin(V, minvotes = 99), "Too few roll calls")
  expect_error(anomin(V, minvotes = 10, polarity = 99), "`polarity` must")
})

test_that("ordirt() recovers ideal points and keeps each item's categories", {
  # the simulation of rmorie's own test of the sampler: 120 respondents,
  # 8 items of 4 categories, 30 answers missing
  set.seed(7)
  th <- rnorm(120)
  eta <- outer(th, runif(8, 0.8, 1.6)) +
    matrix(rnorm(8, 0, 0.5), 120, 8, byrow = TRUE)
  Y <- matrix(cut(eta + rnorm(960), c(-Inf, -0.5, 0.4, 1.2, Inf),
    labels = FALSE), 120, 8)
  Y[sample(length(Y), 30)] <- NA
  f <- ordirt(Y, nsamp = 300, burn = 300, seed = 1)
  expect_gt(cor(f$ideal[, 1], th), 0.9)
  expect_gt(f$disc[1, 1], 0)
  expect_true(all(vapply(f$cuts, function(g) g[1] == 0 && all(diff(g) > 0),
    NA)))
  Y2 <- Y
  Y2[, 1] <- c(10, 20, 30, 40)[Y[, 1]]
  g <- ordirt(Y2, nsamp = 20, burn = 10)
  expect_equal(g$ideal, ordirt(Y, nsamp = 20, burn = 10)$ideal,
               tolerance = testthat_tolerance())
  expect_error(ordirt(cbind(Y, 1)), "two observed categories")
})

test_that("dynirt() recovers paths; -1/1 votes and bad periods", {
  set.seed(7)
  period <- rep(1:3, each = 15)
  walk <- function(s) s + cumsum(c(0, rnorm(2, 0, 0.3)))
  truth <- t(vapply(rnorm(12), walk, numeric(3)))
  p <- pnorm(-rep(rnorm(45, 0, 0.5), each = 12) +
    rep(rnorm(45, 1.5, 0.5), each = 12) * truth[, period])
  V <- matrix(rbinom(length(p), 1, p), 12)
  f <- dynirt(V, period, nsamp = 100, burn = 50)
  expect_gt(abs(cor(as.numeric(f$theta), as.numeric(truth))), 0.8)
  expect_identical(dim(f$theta), c(12L, 3L))
  g <- dynirt(2 * V - 1, period, nsamp = 100, burn = 50)
  expect_equal(g$theta, f$theta, tolerance = testthat_tolerance())
  h <- dynirt(V, as.character(period), nsamp = 10, burn = 5)
  expect_identical(h$periods, c("1", "2", "3"))
  expect_error(dynirt(V, 1:2), "one period per roll call")
  expect_error(dynirt(V * 3, period), "0/1")
})

test_that("cjrirt() from its own EM start is the default chain", {
  set.seed(4)
  v <- matrix(rbinom(80, 1, plogis(outer(rnorm(8), rnorm(10)) * 3)), 8)
  v[2, 3] <- NA
  # the default start is irtem() with its defaults
  a <- cjrirt(v, 1, 10, 5, seed = 3)
  b <- cjrirt(v, 1, 10, 5, seed = 3, start = irtem(v)$ideal)
  expect_identical(b$ideal, a$ideal)
  expect_identical(b$beta, a$beta)
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(cjrirt(v, 1, 10, 5, start = 1:3),
     "cjrirt: `start` must hold 8 ideal points.")
  ef(cjrirt(v, 1, 10, 5, start = "a"),
     "cjrirt: `start` must be a numeric matrix.")
  ef(cjrirt(v, chains = NA), "cjrirt: `chains` must be TRUE or FALSE.")
  # a flat prior on the roll-call parameters and identical ideal points:
  # the intercept and the slope cannot be told apart
  ef(cjrirt(v, 1, 10, 5, bvar = Inf, start = rep(0, 8)),
     paste("cjrirt: the conditional covariance of the roll-call parameters",
          "is not positive definite."))
  # ideal points at +-60 put the votes against them 30 standard
  # deviations into the tail: the truncated draws stay finite
  st <- ifelse(rowMeans(v, na.rm = TRUE) > 0.5, 60, -60)
  st[1] <- -st[1]
  x <- cjrirt(v, 1, 20, 0, start = st, seed = 1)
  expect_true(all(is.finite(x$xchain)))
  expect_true(all(is.finite(x$alpha)) && all(is.finite(x$beta)))
})

test_that("bamscl() puts the polarity stimulus on the left", {
  set.seed(2)
  tr <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
  Z <- t(sapply(1:12, function(i) {
    4 + rnorm(1, 0, 0.3) + runif(1, 0.6, 1.4) * tr + rnorm(5, 0, 0.3)
  }))
  # the fifth stimulus is the rightmost: as the polarity it is flipped
  # to the left, and the order of the others with it
  r <- bamscl(Z, 10, 10, polarity = 5)
  expect_lt(r$zeta[5], 0)
  expect_identical(order(r$zeta), 5:1)
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(bamscl(Z, 1), "bamscl: `nsamp` must be one whole number of at least 2.")
  ef(bamscl(Z, polarity = 6),
     "bamscl: `polarity` must be a column of `Z`, at most 5.")
})

test_that("the dimensions of the samplers are checked", {
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  set.seed(7)
  D <- as.matrix(dist(matrix(rnorm(12), 6)))
  ef(bmdsf(D, 0), "bmdsf: `ndim` must be one whole number of at least 1.")
  ef(bmdsf(D, 6), "bmdsf: `ndim` must be below the number of objects, 6.")
  ef(bunfld(D, 0), "bunfld: `ndim` must be one whole number of at least 1.")
  ef(bunfld(D[, 1:3], 4),
     "bunfld: `ndim` must be at most 3, the smaller dimension of `D`.")
  Y <- matrix(sample(1:3, 60, TRUE), 20)
  ef(ordirt(Y, 0), "ordirt: `ndim` must be one whole number of at least 1.")
  ef(ordirt(Y, 4),
     "ordirt: `ndim` must be at most 3, the smaller dimension of `Y`.")
  # two dimensions: the draws are rotated onto the first, which is kept
  o2 <- ordirt(Y, ndim = 2, nsamp = 5, burn = 2)
  expect_identical(dim(o2$ideal), c(20L, 2L))
  # an item two respondents answered, with no prior on its intercept and
  # two loadings: three coefficients from two answers
  Y[3:20, 1] <- NA
  Y[1:2, 1] <- 1:2
  ef(ordirt(Y, ndim = 2, nsamp = 5, burn = 2),
     paste("ordirt: the conditional covariance of an item's loadings is not",
           "positive definite."))
})

test_that("anomin() refuses polarities and too few roll calls", {
  set.seed(4)
  x <- seq(-1, 1, length.out = 20)
  V <- outer(x, runif(40, -0.7, 0.7), ">") * 1
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(anomin(V, polarity = "a", minvotes = 10),
     "anomin: `polarity` must be a numeric vector.")
  ef(anomin(V, ndim = 2, polarity = 1:3, minvotes = 10),
     "anomin: `polarity` must give one legislator per dimension, 2.")
  # two roll calls cannot place legislators in three dimensions
  ef(anomin(V[, c(10, 30)], ndim = 3, minvotes = 1, polarity = 20),
     "anomin: the cross-product of the coordinates is singular.")
})

test_that("dynirt() anchors, and refuses what it cannot fit", {
  set.seed(7)
  period <- rep(1:3, each = 15)
  walk <- function(s) s + cumsum(c(0, rnorm(2, 0, 0.3)))
  truth <- t(sapply(rnorm(12), walk))
  p <- pnorm(-rep(rnorm(45, 0, 0.5), each = 12) +
    rep(rnorm(45, 1.5, 0.5), each = 12) * truth[, period])
  V <- matrix(rbinom(length(p), 1, p), 12)
  # anchored at the legislator furthest to the left: the anchor's mean is
  # made positive, so the paths come out reversed
  lo <- which.min(rowMeans(truth))
  f <- dynirt(V, period, nsamp = 100, burn = 50, anchor = lo)
  expect_gt(mean(f$theta[lo, ]), 0)
  expect_lt(cor(as.numeric(f$theta), as.numeric(truth)), -0.8)
  ef <- function(call, msg) expect_error(call, msg, fixed = TRUE)
  ef(dynirt(V, c(1, NA, period[-(1:2)]), nsamp = 5, burn = 2),
     "dynirt: `period` must not be NA.")
  ef(dynirt(V, period, nsamp = 5, burn = 2, anchor = 13),
     "dynirt: `anchor` must be a row of `votes`, at most 12.")
  ef(dynirt(V, period, nsamp = 5, burn = 2, anchor = "a"),
     "dynirt: `anchor` must be one whole number of at least 1.")
  # one legislator and flat priors: an intercept and a slope from one
  # ideal point per period
  ef(dynirt(V[1, , drop = FALSE], period, nsamp = 5, burn = 2, A0 = 0,
            B0 = 0),
     paste("dynirt: the conditional covariance of a roll call's parameters",
           "is not positive definite."))
  # every legislator votes alike: no first principal component to start
  # from, so every start is 0
  same <- matrix(V[1, ], 12, 45, byrow = TRUE)
  g <- dynirt(same, period, nsamp = 5, burn = 2)
  expect_identical(dim(g$theta), c(12L, 3L))
  expect_true(all(is.finite(g$theta)))
})
