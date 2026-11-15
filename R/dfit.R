#' Fit and rank candidate distributions
#'
#' Fits candidate distributions to one numeric variable by maximum
#' likelihood, ranks them by AIC, and reports for each the BIC, a
#' Kolmogorov-Smirnov distance, an Anderson-Darling statistic and percentile
#' bootstrap intervals for its parameters. A Hill estimate of the tail index
#' is added for positive data. Use it to choose the model of a count of
#' offences, a sentence length or a time to reoffending before modelling it.
#'
#' @details The families and their parameters: \code{normal} (\code{mean},
#'   \code{sd}), \code{lognormal} (\code{meanlog}, \code{sdlog}),
#'   \code{exponential} (\code{rate}), \code{gamma} (\code{shape},
#'   \code{rate}), \code{weibull} (\code{shape}, \code{scale}), \code{t}
#'   (\code{location}, \code{scale}, \code{df}: a location-scale Student t),
#'   \code{pareto} (\code{xmin}, \code{alpha}, with \code{xmin} the smallest
#'   value, its maximum likelihood estimate), \code{beta} (\code{shape1},
#'   \code{shape2}) for proportions, and for counts \code{poisson}
#'   (\code{lambda}), \code{negbin} (\code{size}, \code{mu}), \code{zip}
#'   (\code{pi}, \code{lambda}) and \code{zinb} (\code{pi}, \code{size},
#'   \code{mu}), the zero-inflated forms with \code{pi} the share of
#'   structural zeros. The normal, lognormal, exponential, Pareto and Poisson
#'   estimates are closed form; the negative binomial searches its size alone
#'   (its mu estimate is the mean) by Brent's method on the log scale over
#'   (-10, 15), as \code{optimize()}; the rest are maximised on a log or logit
#'   scale by the Nelder-Mead simplex (at most 4000 evaluations, relative
#'   tolerance 1e-12) and then BFGS (at most 500 iterations, relative
#'   tolerance 1e-14, central-difference gradients), the algorithms of
#'   \code{optim()}, in the compiled core.
#'
#' With \code{fams = NULL} every family the data allow is fitted: for whole
#' numbers of 0 or more the count families; otherwise the positive families
#' for values above 0 (\code{exponential} for values of 0 or more) and
#' \code{beta} for values strictly between 0 and 1; and \code{normal} and
#' \code{t} always. A family named in \code{fams} that the data do not allow
#' is refused. A fit that fails (the optimiser stops at a non-finite
#' likelihood) is a row with \code{converged = FALSE} and no statistics.
#'
#' AIC is \eqn{2k - 2\ell} and BIC \eqn{k \log n - 2\ell}. The
#' Kolmogorov-Smirnov and Anderson-Darling statistics use the fitted
#' parameters, so their usual null distributions do not apply; they rank
#' fits, they are not tests. The bootstrap resamples the data \code{nboot}
#' times and refits; the interval is the percentile interval at
#' \code{conflvl}. The Hill estimator uses the \code{k} largest values
#' (default \eqn{\lfloor\sqrt n\rfloor}): \eqn{\hat\alpha = k / \sum_{i \le
#' k} \log(x_{(n-i+1)} / x_{(n-k)})}, with the asymptotic interval
#' \eqn{\hat\alpha (1 \pm z / \sqrt k)}.
#'
#' Missing values are dropped and counted (\code{nmiss}). Errors: \code{x}
#' not numeric, an infinite value, fewer than five values, all values equal,
#' an unknown family or one the data do not allow. The bootstrap is drawn
#' from R's generator seeded with \code{seed} for this call only; the
#' session's random stream is left as it was.
#'
#' @param x Numeric vector, the variable.
#' @param fams Character vector of families (see Details), or \code{NULL}
#'   for every family the data allow.
#' @param nboot Whole number, default 200: bootstrap resamples (0 for none).
#' @param conflvl Number in (0, 1), default 0.95: the bootstrap and Hill
#'   interval level.
#' @param seed Whole number, default 1: the seed of the bootstrap.
#' @param k Whole number or \code{NULL}: the upper order statistics of the
#'   Hill estimate.
#' @param ... Not used.
#' @return A list of class \code{c("dstfit", "list")}:
#'   \describe{
#'     \item{\code{table}}{data frame, one row per family ranked by AIC:
#'       \code{rank}, \code{family}, \code{npar}, \code{loglik}, \code{aic},
#'       \code{bic}, \code{ks}, \code{ad}, \code{converged} and
#'       \code{params} (the estimates as text).}
#'     \item{\code{fits}}{named list of the estimates (named numeric vectors)
#'       of the converged fits.}
#'     \item{\code{boot}}{data frame: \code{family}, \code{param},
#'       \code{estimate}, \code{lower}, \code{upper} (\code{NA} without a
#'       bootstrap).}
#'     \item{\code{hill}}{list: \code{k}, \code{alpha}, \code{lower},
#'       \code{upper}, or \code{NULL} for data that are not all positive or
#'       fewer than ten values.}
#'     \item{\code{x}, \code{nmiss}, \code{discrete}, \code{conflvl}}{the
#'       data used (sorted), the count dropped, whether the data are counts,
#'       and the level.}
#'   }
#' @references Hill, B. M. (1975). A simple general approach to inference
#'   about the tail of a distribution. Annals of Statistics 3, 1163-1174.
#'   \doi{10.1214/aos/1176343247}
#'
#'   Anderson, T. W. and Darling, D. A. (1954). A test of goodness of fit.
#'   Journal of the American Statistical Association 49, 765-769.
#'   \doi{10.1080/01621459.1954.10501232}
#'
#'   Nelder, J. A. and Mead, R. (1965). A simplex method for function
#'   minimization. The Computer Journal 7(4), 308-313.
#'   \doi{10.1093/comjnl/7.4.308}
#' @seealso \code{\link{profcols}} for the column's summaries.
#' @examples
#' set.seed(1)
#' x <- stats::rgamma(300, shape = 2, rate = 0.5) # simulated
#' f <- dstfit(x, nboot = 20)
#' f
#' f$fits$gamma
#'
#' # simulated counts with extra zeros
#' n <- ifelse(stats::runif(300) < 0.3, 0, stats::rpois(300, 4))
#' dstfit(n, nboot = 0)$table[, c("family", "aic")]
#'
#' # a family the data do not allow
#' try(dstfit(c(-1, 2, 3, 4, 5), fams = "gamma"))
#' @export
dstfit <- function(x, fams = NULL, nboot = 200L, conflvl = 0.95, seed = 1L,
                   k = NULL) {
  f <- "dstfit"
  if (!is.numeric(x) || is.object(x)) {
    .ckstop(f, "`x` must be a numeric vector")
  }
  if (any(is.infinite(x))) .ckstop(f, "`x` has an infinite value")
  nmiss <- sum(is.na(x))
  x <- sort(as.numeric(x[!is.na(x)]))
  if (length(x) < 5L) {
    .ckstop(f, "`x` needs at least five values that are not missing")
  }
  if (length(unique(x)) < 2L) {
    .ckstop(
      f,
      "every value of `x` is the same; there is no distribution to fit"
    )
  }
  nboot <- .ckn(nboot, "nboot", f, integer = TRUE)
  if (length(nboot) != 1L || is.na(nboot) || nboot < 0 || nboot > 1e5) {
    .ckstop(f, "`nboot` must be one whole number from 0 to 1e+05")
  }
  conflvl <- .ckn(conflvl, "conflvl", f)
  if (length(conflvl) != 1L || is.na(conflvl) || conflvl <= 0 ||
        conflvl >= 1) {
    .ckstop(f, "`conflvl` must be a single number strictly between 0 and 1")
  }
  fm <- .dfams()
  pos <- all(x > 0)
  cnt <- all(x >= 0 & x == round(x))
  okf <- c(TRUE, pos, all(x >= 0), pos, pos, TRUE, pos, all(x > 0 & x < 1),
           cnt, cnt, cnt && any(x == 0), cnt && any(x == 0))
  disc <- cnt
  if (is.null(fams)) {
    fams <- fm$name[okf & (fm$disc | !disc | fm$name %in% c("normal", "t"))]
  }
  fams <- .ckargs(fams, fm$name, "fams", f)
  for (nm in fams) {
    if (!okf[fm$name == nm]) {
      .ckstop(f, "`%s` needs %s", nm, fm$need[fm$name == nm])
    }
  }
  n <- length(x)
  fits <- list()
  rows <- lapply(fams, function(nm) {
    i <- match(nm, fm$name)
    p <- .Call(C_clab_r_dfit, i - 1L, x)
    if (is.null(p)) {
      return(data.frame(family = nm, npar = NA_integer_, loglik = NA_real_,
                        aic = NA_real_, bic = NA_real_, ks = NA_real_,
                        ad = NA_real_, converged = FALSE, params = "",
                        stringsAsFactors = FALSE))
    }
    names(p) <- fm$par[[i]]
    fits[[nm]] <<- p
    g <- .Call(C_clab_r_dfgof, i - 1L, x, p)
    data.frame(
      family = nm, npar = length(p), loglik = g[1L],
      aic = 2 * length(p) - 2 * g[1L], bic = log(n) * length(p) - 2 * g[1L],
      ks = g[2L], ad = g[3L], converged = TRUE,
      params = toString(sprintf("%s = %s", names(p), format(signif(p, 4)))),
      stringsAsFactors = FALSE
    )
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(tab$aic, na.last = TRUE), , drop = FALSE]
  tab <- cbind(rank = seq_len(nrow(tab)), tab)
  rownames(tab) <- NULL
  fits <- fits[intersect(tab$family, names(fits))]
  a <- (1 - conflvl) / 2
  idx <- if (nboot > 0L) {
    seed <- .ckseed(seed, "seed", f)
    .clseed(seed, lapply(names(fits), function(nm) {
      vapply(seq_len(nboot), function(b) sample.int(n, n, replace = TRUE),
             integer(n))
    }))
  }
  boot <- do.call(rbind, c(
    list(data.frame(family = character(), param = character(),
                    estimate = numeric(), lower = numeric(), upper = numeric(),
                    stringsAsFactors = FALSE)),
    lapply(seq_along(fits), function(j) {
      p <- fits[[j]]
      lim <- if (nboot > 0L) {
        .Call(C_clab_r_dfboot, match(names(fits)[j], fm$name) - 1L, x,
              matrix(idx[[j]], nrow = n), a)
      } else {
        matrix(NA_real_, 2L, length(p))
      }
      data.frame(family = names(fits)[j], param = names(p),
                 estimate = unname(p), lower = lim[1L, ], upper = lim[2L, ],
                 stringsAsFactors = FALSE)
    })
  ))
  hill <- NULL
  if (n >= 10L && all(x > 0)) {
    kk <- if (is.null(k)) {
      floor(sqrt(n))
    } else {
      v <- .ckn(k, "k", f, integer = TRUE)
      if (length(v) != 1L || is.na(v) || v < 2 || v > n - 1) {
        .ckstop(f, "`k` must be one whole number from 2 to %s", format(n - 1))
      }
      v
    }
    h <- .Call(C_clab_r_dfhill, x, as.integer(kk), conflvl)
    hill <- list(k = as.integer(kk), alpha = h[1L], lower = h[2L],
                 upper = h[3L])
  }
  structure(
    list(table = tab, fits = fits, boot = boot, hill = hill, x = x,
         nmiss = nmiss, discrete = disc, conflvl = conflvl),
    class = c("dstfit", "list")
  )
}

#' @rdname dstfit
#' @export
print.dstfit <- function(x, ...) {
  t <- x$table
  cat("Distributions fitted by maximum likelihood:", length(x$x), "values\n")
  print(t[c("rank", "family", "npar", "loglik", "aic", "bic", "ks", "ad",
            "params")], row.names = FALSE, digits = 4)
  if (!is.null(x$hill)) {
    cat(sprintf("Hill tail index (k = %d): %.3g (%.3g to %.3g)\n", x$hill$k,
                x$hill$alpha, x$hill$lower, x$hill$upper))
  }
  invisible(x)
}

# The families: names, whether they are for counts, what they need, their
# parameters.
.dfams <- function() {
  pos <- "every value above 0"
  cn <- "counts (whole numbers, 0 or more)"
  z <- "counts with at least one zero"
  list(
    name = c("normal", "lognormal", "exponential", "gamma", "weibull", "t",
             "pareto", "beta", "poisson", "negbin", "zip", "zinb"),
    disc = c(rep(FALSE, 8L), rep(TRUE, 4L)),
    need = c("numbers", pos, "no value below 0", pos, pos, "numbers", pos,
             "every value between 0 and 1, exclusive", cn, cn, z, z),
    par = list(
      c("mean", "sd"), c("meanlog", "sdlog"), "rate", c("shape", "rate"),
      c("shape", "scale"), c("location", "scale", "df"), c("xmin", "alpha"),
      c("shape1", "shape2"), "lambda", c("size", "mu"), c("pi", "lambda"),
      c("pi", "size", "mu")
    )
  )
}
