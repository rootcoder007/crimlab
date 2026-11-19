#' ppmlcf(family = "nb"): the control function with a negative binomial
#' second stage (Feigenbaum and Muller 2016, after Cameron and Trivedi
#' 2013): the fixed effects as indicator columns (the first level of each
#' left out, an intercept kept), the first stage by least squares of each
#' endogenous column on the instruments, the controls and the indicators,
#' its residuals added to the second stage, fitted by maximum likelihood
#' with theta (clab_cntreg()); the errors by a bootstrap of the clusters
#' refitting both stages, or naive (the clustered sandwich of the second
#' stage alone).
#' @noRd
.axcfnb <- function(d, yv, endo, inst, x, fe, off, clus, cl, se, ssc, nboot,
                    seed, level, maxit, f) {
  tm <- c(endo, x, paste0("vhat:", endo))
  p <- length(tm)
  fit <- function(rows, cc = NULL, ncc = 0L) {
    dd <- d[rows, , drop = FALSE]
    dm <- .axcfdum(dd, fe)
    zz <- cbind(1, .axmat(dd, c(inst, x), f), dm)
    vh <- vapply(endo, function(e) {
      stats::lm.fit(zz, as.numeric(dd[[e]]))$residuals
    }, numeric(nrow(dd)))
    x2 <- cbind(.axmat(dd, c(endo, x), f), matrix(vh, nrow(dd)), 1, dm)
    .Call(C_clab_r_nbfit, as.numeric(yv[rows]), x2,
          if (is.null(off)) NULL else off[rows], cc, as.integer(ncc),
          as.integer(maxit), 1e-10)
  }
  if (se == "boot") {
    r <- fit(seq_along(yv))
    b <- .axboot(d, clus[1L], nboot, seed, p, function(rows, w) {
      fit(rep(rows, w))[[1L]][seq_len(p)]
    }, f)
    sev <- b$se
    dr <- b$draws
    nf <- b$nfail
  } else {
    r <- fit(seq_along(yv), cl$c1, cl$n1)
    g <- cl$n1
    n <- length(yv)
    k <- ncol(r[[2L]])
    fac <- switch(ssc, G = g / (g - 1),
                  nk = g / (g - 1) * (n - 1) / (n - k), none = 1)
    sev <- sqrt(fac * diag(r[[3L]])[seq_len(p)])
    dr <- NULL
    nf <- 0L
  }
  tab <- .axtab(tm, r[[1L]][seq_len(p)], sev, level)
  tab$rr <- exp(tab$est)
  list(coef = tab, draws = dr, theta = r[[4L]][3L], nfail = nf)
}

#' The indicators of the fixed effects of the rows of d, the first level
#' of each left out.
#' @noRd
.axcfdum <- function(d, fe) {
  m <- lapply(fe, function(v) {
    stats::model.matrix(~g, data.frame(g = factor(d[[v]])))[, -1L,
                                                           drop = FALSE]
  })
  do.call(cbind, c(list(matrix(0, nrow(d), 0)), m))
}
