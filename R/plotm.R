#' Plot methods for the package's results
#'
#' Every result of the package that holds something to draw has a `plot()`
#' method, on base graphics only. No argument beyond the object is needed;
#' every method takes the same customisation: `main`, `col` or `pal`, and
#' `...` passed to the base-graphics call, with the session look set by
#' [plotopt()]. Each method returns, invisibly, the data it drew (a data
#' frame, or a list of them), so the same picture can be drawn with any
#' other tool. Series that differ in colour also differ in line type or
#' point shape, and shaded cells carry their numbers, so no plot depends on
#' colour alone.
#'
#' @details What each method draws:
#' * Estimates with intervals, as a dot-and-interval (forest) chart:
#'   `cntreg` (rate ratios, log scale), `dastd`, `grprr`, `sirex`, `svymx`,
#'   `mndla`, `thrtst` (thresholds), `sarml` and `spcure` (coefficients),
#'   `mtrbnd` (bounds by treatment), `pbbnd` (bounds of published cells),
#'   `dkfig` (bounds on the dark figure), `didfe` (with its Goodman-Bacon
#'   comparisons and leave-one-out estimates when computed), `smdbal`
#'   (covariate balance), `jmono` (rates by judge), `relchk` (the readings
#'   of a table against the claim) and `dissix` (segregation indices).
#' * Effects by time since adoption: `evtst` and `stgdid`, before adoption
#'   in open circles and after in filled ones.
#' * Paths over time: `synthc` and `sdidx` (treated against its synthetic
#'   comparison), `trndrr` (series and fitted trend), `stkflw` (population,
#'   person-days, average daily population and length of stay), `dktflw`
#'   (clearance rate and pending cases), `jpipe`, `hsper`, `yoypx` (percent
#'   change with its interval), `cumrsk`, `agemix` and `apcie`.
#' * Curves: `kmrec` and `lftab` (survival), `lnroc` (lineup ROC curves),
#'   `lnkpr` (the ROC of the linkage model), `cconc` (Lorenz curve),
#'   `crossk` (cross K function), `bandx` (the statistic against the top
#'   band's cap), `jdgiv` (sensitivity of the instrument), `xcser`
#'   (cross-correlations) and `hwkres` (time-rescaling residuals).
#' * Bars: `arwt`, `crsev`, `gmixlr` and `pedlr` (log10 likelihood ratio
#'   by locus), `hypq` (workload of each unit), `lnsize`, `lostay`,
#'   `mmshr` (stacked, hatched), `pinpl` and `rptvic` (observed against
#'   expected).
#' * Matrices and surfaces: `ctabx` (residuals), `knxnr` (Knox ratios),
#'   `qdrtst`, `trnmat`, `zcomp`, `kdsurf`, `prohs` (hot cells outlined)
#'   and `sdbsim`.
#' * Points: `bufas` (distances to sites), `cmpsc` (segment positions),
#'   `ebshr` (raw against shrunk values), `ibench` (rate against
#'   benchmark) and `kscan` (clusters as circles).
#'
#' @param x A result of the package.
#' @inheritParams plotlrnz
#' @return Invisibly, the data drawn: a data frame, or for a method of
#'   several panels a list of data frames. A result missing a part the plot
#'   needs is refused.
#' @srrstats {EA5.0a} Text follows plotopt(cex = ); no typeface is set.
#' @srrstats {EA5.4} Axes are drawn by base graphics from pretty().
#' @name plotm
#' @examples
#' d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
#' i <- match(d$town, LETTERS)
#' first <- rep(c(2019, 2021, NA), each = 4)
#' d$first <- first[i]
#' d$bylaw <- as.integer(!is.na(d$first) & d$year >= d$first)
#' d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * d$bylaw + sin(seq_len(120))
#' plot(evtst(d$rate, d$town, d$year, d$first, leads = 3, lags = 3))
#' plot(didfe(d$rate, d$town, d$year, d$bylaw, bacon = TRUE, loo = TRUE))
#' one <- d[d$town %in% c("A", "I", "J", "K", "L"), ]
#' plot(synthc(one$rate, one$town, one$year, one$first))
#' plot(grprr(c(a = 120, b = 95, c = 30), c(a = 6e4, b = 2.1e4, c = 1.5e4)))
#' plot(kmrec(c(3, 8, 8, 12, 20, 24), event = c(1, 1, 0, 1, 0, 1),
#'            grp = c("a", "a", "a", "b", "b", "b")))
#' plot(trnmat(c("N", "N", "E", "N", "W", "W", "E", "E"),
#'             c(1, 1, 1, 2, 2, 2, 3, 3), time = c(1, 2, 3, 1, 2, 3, 1, 2)))
#' set.seed(1)
#' plot(cconc(rpois(200, 2), nsim = 0))
NULL

# ---- estimates with intervals ------------------------------------------------

#' @rdname plotm
#' @export
plot.cntreg <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cntreg", "irr")
  r <- x$irr
  .pdots(r$term, r$irr, r$lo, r$hi,
         .pmain(main, "Rate ratios with their intervals"),
         "Rate ratio (log scale)", col, pal, ref = 1,
         grp = if (length(unique(r$fit)) > 1L) r$fit, log = "x", ...)
}

#' @rdname plotm
#' @export
plot.dastd <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dastd", "rates")
  r <- x$rates
  .pdots(r$group, r$rate, r$lo, r$hi,
         .pmain(main, "Directly standardised rates"),
         paste("Rate per", format(x$per)), col, pal, ...)
}

#' @rdname plotm
#' @export
plot.grprr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "grprr", "rates")
  r <- x$rates
  .pdots(r$group, r$rate, r$lo, r$hi,
         .pmain(main, "Rates by group with their intervals"),
         if (identical(x$type, "rate")) paste("Rate per", format(x$per))
         else "Share", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.sirex <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "sirex", "sir")
  r <- x$sir
  .pdots(r$group, r$sir, r$lo, r$hi,
         .pmain(main, "Standardised ratios with their intervals"),
         "Observed / expected", col, pal, ref = 1,
         grp = ifelse(r$out & !is.na(r$out), "outside the funnel limits",
                      "inside the funnel limits"), ...)
}

#' @rdname plotm
#' @export
plot.svymx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "svymx", "est")
  r <- x$est
  .pdots(r$dom, r$est, r$lo, r$hi,
         .pmain(main, "Design-based estimates with their intervals"),
         if (is.null(x$stat)) "Estimate" else x$stat, col, pal, ...)
}

#' @rdname plotm
#' @export
plot.mndla <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "mndla", "tab")
  r <- x$tab
  .pdots(r$group, r$share, r$lo, r$hi,
         .pmain(main, "Share of placements past the limit"), "Share", col,
         pal, ...)
}

#' @rdname plotm
#' @export
plot.thrtst <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "thrtst", "thr")
  r <- x$thr
  .pdots(r$group, r$thr, r$lo, r$hi,
         .pmain(main, "Search thresholds by group"), "Threshold", col, pal,
         ...)
}

#' @rdname plotm
#' @export
plot.sarml <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "sarml", "coef")
  r <- x$coef
  z <- stats::qnorm(0.975)
  .pdots(r$term, r$est, r$est - z * r$se, r$est + z * r$se,
         .pmain(main, "Coefficients with 95% intervals"), "Estimate", col,
         pal, ref = 0, ...)
}

#' @rdname plotm
#' @export
plot.spcure <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "spcure", "coef")
  r <- x$coef
  .pdots(paste0(r$part, ": ", r$term), r$est, r$lo, r$hi,
         .pmain(main, "Cure model coefficients with their intervals"),
         "Estimate", col, pal, ref = 0, grp = r$part, ...)
}

#' @rdname plotm
#' @export
plot.mtrbnd <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "mtrbnd", "bnd")
  r <- x$bnd
  .pdots(as.character(r$z), r$mean, r$lower, r$upper,
         .pmain(main, "Observed means (points) and bounds by treatment"),
         "Mean outcome", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.pbbnd <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "pbbnd", "bnd")
  r <- x$bnd
  .pdots(sprintf("cell %d", seq_len(nrow(r))), r$value, r$lower, r$upper,
         .pmain(main, "Published values (points) and their bounds"),
         "Value", col, pal, grp = r$status, ...)
}

#' @rdname plotm
#' @export
plot.dkfig <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dkfig", "bnd")
  r <- x$bnd
  .pdots(sprintf("v %g, r %g", r$v, r$r), r$v, r$lo, r$hi,
         .pmain(main, "Recorded rates (points) and bounds on the true rate"),
         "Rate", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.smdbal <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "smdbal", "tab")
  r <- x$tab[order(-abs(x$tab$smd)), , drop = FALSE]
  d <- .pdots(r$var, abs(r$smd), NA, NA, .pmain(main, "Covariate balance"),
              "Absolute standardised difference", col, pal, ref = x$thr,
              ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.jmono <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "jmono", "judges")
  r <- x$judges[order(-x$judges$rate), , drop = FALSE]
  .pdots(r$judge, r$rate, NA, NA, .pmain(main, "Treatment rate by judge"),
         "Rate", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.relchk <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "relchk", "reads")
  r <- x$reads
  .pdots(sprintf("%s %s: %s vs %s, %s", r$grp, r$type, r$group, r$ref,
                 r$event),
         r$value, NA, NA,
         .pmain(main, sprintf("Readings of the table; the claim is %s",
                              format(x$value))),
         "Value", col, pal, ref = x$value,
         grp = ifelse(r$match, "reproduces the claim", "does not"), ...)
}

#' @rdname plotm
#' @export
plot.dissix <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dissix", "grp")
  r <- x$grp
  k <- c(is = "isolation", h = "entropy", atk = "Atkinson", gini = "Gini")
  .pdots(rep(r$group, 4L), unlist(r[names(k)], use.names = FALSE), NA, NA,
         .pmain(main, "Segregation indices by group"), "Index", col, pal,
         grp = rep(k, each = nrow(r)), ...)
}

#' @rdname plotm
#' @export
plot.didfe <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "didfe", "coef")
  n <- 1L + sum(!vapply(x[c("bacon", "loo")], is.null, NA))
  op <- .ppanels(n)
  on.exit(graphics::par(op), add = TRUE)
  r <- x$coef
  out <- list(coef = .pdots(r$term, r$est, r$lo, r$hi,
                            .pmain(main, "Two-way fixed-effects estimate"),
                            "Effect", col, pal, ref = 0, ...))
  if (!is.null(x$bacon)) {
    b <- x$bacon
    g <- factor(b$type, levels = unique(b$type))
    cl <- .pcols(nlevels(g), col, pal)
    .ppar()
    .pfrm(.prange(b$weight, 0), .prange(b$est, r$est[1L]),
            "Goodman-Bacon comparisons", "Weight", "Estimate", ...)
    graphics::abline(h = r$est[1L], lty = 2, col = "#636363")
    graphics::points(b$weight, b$est, col = cl[g], pch = .ppch(g), cex = 1.3)
    .pkey("bottomright", levels(g), cl, pch = .ppch(seq_len(nlevels(g))))
    out$bacon <- b
  }
  if (!is.null(x$loo)) {
    l <- x$loo
    z <- stats::qnorm(1 - (1 - x$level) / 2)
    out$loo <- .pdots(paste("without", l$unit), l$est, l$est - z * l$se,
                      l$est + z * l$se, "Leaving one adopting unit out",
                      "Effect", col, pal, ref = r$est[1L], ...)
  }
  invisible(out)
}

# ---- effects by time since adoption -----------------------------------------

#' @rdname plotm
#' @export
plot.evtst <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "evtst", "est")
  e <- .pev(x)
  .ppath(e$evtm, e$est, e$lo, e$hi, .pmain(main, "Event study"),
         "Periods since adoption", "Effect, with its interval", col, pal,
         ...)
}

#' @rdname plotm
#' @export
plot.stgdid <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "stgdid", "dynamic")
  e <- .pev(x)
  .ppath(e$evtm, e$est, e$lo, e$hi,
         .pmain(main, "Effect by time since adoption"),
         "Periods since adoption", "Effect, with its interval", col, pal,
         ...)
}

# ---- paths over time --------------------------------------------------------

#' @rdname plotm
#' @export
plot.synthc <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "synthc", "path")
  p <- x$path
  op <- .ppanels(2L)
  on.exit(graphics::par(op), add = TRUE)
  .plines(p$time, list(treated = p$treated, synthetic = p$synth),
          .pmain(main, paste("Synthetic control:", x$unit)), "Period",
          "Outcome", col, pal, where = "topleft", ...)
  .plines(p$time, list(gap = p$gap), "Treated minus synthetic", "Period",
          "Gap", col, pal, ref = 0, ...)
  invisible(p)
}

#' @rdname plotm
#' @export
plot.sdidx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "sdidx", "path")
  p <- x$path
  .plines(p$time, list(treated = p$treated, `weighted comparison` = p$synth),
          .pmain(main, "Synthetic difference in differences"), "Period",
          "Outcome", col, pal, where = "topleft", ...)
  invisible(p)
}

#' @rdname plotm
#' @export
plot.trndrr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "trndrr", c("yoy", "fitted", "rr"))
  d <- data.frame(x = x$yoy$x, y = x$yoy$y, fitted = x$fitted)
  .plines(d$x, list(observed = d$y, fitted = d$fitted),
          .pmain(main, sprintf("Trend: rate ratio %.3f per period",
                               x$rr$rr[1L])), "Period", "Value", col, pal,
          ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.stkflw <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "stkflw", "flow")
  f <- x$flow
  k <- c(people = "People held", days = "Person-days",
         adp = "Average daily population", alos = "Average length of stay")
  main <- .pmain(main, "")
  op <- .ppanels(4L)
  on.exit(graphics::par(op), add = TRUE)
  cl <- .pcols(4L, col, pal)
  for (i in seq_along(k)) {
    .plines(f$period, stats::setNames(list(f[[names(k)[i]]]), names(k)[i]),
            k[[i]], "Period", "", cl[i], ...)
  }
  if (nzchar(main)) graphics::mtext(main, side = 3, line = -1.2, outer = TRUE)
  invisible(f[, c("period", names(k))])
}

#' @rdname plotm
#' @export
plot.dktflw <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dktflw", "flow")
  f <- x$flow
  op <- .ppanels(2L)
  on.exit(graphics::par(op), add = TRUE)
  cl <- .pcols(2L, col, pal)
  .plines(f$period, list(clearance = f$clear),
          .pmain(main, "Clearance rate (disposed / received)"), "Period",
          "Clearance rate", cl[1L], lo = list(f$lo), hi = list(f$hi),
          ref = 1, ...)
  .plines(f$period, list(pending = f$pend), "Pending cases", "Period",
          "Cases", cl[2L], ...)
  invisible(f[, c("period", "clear", "lo", "hi", "pend")])
}

#' @rdname plotm
#' @export
plot.jpipe <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "jpipe", c("index", "rate"))
  d <- data.frame(group = x$group, index = x$index, rate = x$rate)
  .pgrp(d$index, d$rate, d$group,
        .pmain(main, "Rate between stages over time"), "Period", "Rate",
        col, pal, ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.hsper <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "hsper", "periods")
  p <- x$periods
  .plines(p$period, list(`share of events in hot cells` = p$hotshr,
                         `hot cells kept (Jaccard)` = p$jacc),
          .pmain(main, "Hot spots over the periods"), "Period", "Share",
          col, pal, ...)
  invisible(p[, c("period", "hotshr", "jacc")])
}

#' @rdname plotm
#' @export
plot.yoypx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "yoypx", c("prd", "pctchg"))
  d <- data.frame(group = x$group, prd = x$prd, pct = x$pctchg,
                  lo = x$pctlo, hi = x$pctup, stringsAsFactors = FALSE)
  d <- d[is.finite(d$pct), , drop = FALSE]
  g <- factor(d$group, levels = unique(d$group))
  k <- nlevels(g)
  cl <- .pcols(k, col, pal)
  at <- d$prd + (as.integer(g) - (k + 1) / 2) * 0.6 / max(k, 1)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(at), .prange(d$pct, d$lo, d$hi, 0),
          .pmain(main, "Percent change from the period before"), "Period",
          "Percent change", ...)
  graphics::abline(h = 0, lty = 2, col = "#636363")
  w <- is.finite(d$lo) & is.finite(d$hi)
  graphics::segments(at[w], d$lo[w], at[w], d$hi[w], col = cl[g][w], lwd = 2)
  graphics::points(at, d$pct, col = cl[g], pch = .ppch(g), cex = 1.2)
  if (k > 1L) .pkey("topright", levels(g), cl, pch = .ppch(seq_len(k)))
  invisible(d)
}

#' @rdname plotm
#' @export
plot.cumrsk <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cumrsk", c("age", "cumrisk"))
  d <- data.frame(age = x$age, cumrisk = x$cumrisk, lo = x$lo, hi = x$hi)
  .plines(d$age, list(cumrisk = d$cumrisk),
          .pmain(main, "Cumulative risk by age"), "Age", "Cumulative risk",
          col, pal, lo = list(d$lo), hi = list(d$hi), ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.agemix <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "agemix", "curve")
  d <- x$curve
  .plines(d$age, d[setdiff(names(d), "age")],
          .pmain(main, "Age-crime curve and its offender types"), "Age",
          "Rate", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.apcie <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "apcie", "coef")
  d <- x$coef[x$coef$term != "(Intercept)", , drop = FALSE]
  tm <- unique(d$term)
  main <- .pmain(main, "")
  op <- .ppanels(length(tm))
  on.exit(graphics::par(op), add = TRUE)
  for (k in tm) {
    s <- d[d$term == k, ]
    p <- suppressWarnings(as.numeric(s$level))
    num <- !anyNA(p)
    if (!num) p <- seq_len(nrow(s))
    .plines(p, list(est = s$est), paste("Effect of", k), k, "Effect", col,
            pal, lo = list(s$lo), hi = list(s$hi), ref = 0,
            xaxt = if (num) "s" else "n", ...)
    if (!num) graphics::axis(1, at = p, labels = s$level)
  }
  if (nzchar(main)) graphics::mtext(main, side = 3, line = -1.2, outer = TRUE)
  invisible(d)
}

# ---- curves -----------------------------------------------------------------

#' @rdname plotm
#' @export
plot.kmrec <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "kmrec", "km")
  k <- x$km
  s <- unique(k$grp)
  d <- rbind(data.frame(grp = s, time = 0, surv = 1, lo = 1, hi = 1,
                        stringsAsFactors = FALSE),
             k[, c("grp", "time", "surv", "lo", "hi")])
  d <- d[order(match(d$grp, s), d$time), ]
  .pgrp(d$time, d$surv, d$grp, .pmain(main, "Kaplan-Meier survival"),
        "Time", "Share without the event", col, pal, lo = d$lo, hi = d$hi,
        type = "s", ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.lftab <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "lftab", "tab")
  t <- x$tab[is.finite(x$tab$upper), , drop = FALSE]
  d <- data.frame(time = c(x$tab$lower[1L], t$upper), surv = c(1, t$surv),
                  lo = c(1, t$lo), hi = c(1, t$hi))
  .plines(d$time, list(survival = d$surv),
          .pmain(main, "Life table survival"), "Time", "Share surviving",
          col, pal, lo = list(d$lo), hi = list(d$hi), type = "s", ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.lnroc <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "lnroc", "roc")
  r <- x$roc
  s <- unique(r$proc)
  d <- rbind(data.frame(proc = s, far = 0, hr = 0, stringsAsFactors = FALSE),
             r[, c("proc", "far", "hr")])
  d <- d[order(match(d$proc, s), d$far, d$hr), ]
  .pgrp(d$far, d$hr, d$proc, .pmain(main, "Lineup ROC curves"),
        "False identification rate", "Correct identification rate", col,
        pal, where = "bottomright", ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.lnkpr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "lnkpr", "pairs")
  p <- x$pairs
  l <- as.logical(p$linked)
  if (!any(l) || all(l)) {
    stop("plot: an ROC curve needs both linked and unlinked pairs.",
         call. = FALSE)
  }
  th <- sort(unique(p$prob), decreasing = TRUE)
  d <- data.frame(fpr = c(0, vapply(th, function(v) mean(p$prob[!l] >= v),
                                    1)),
                  tpr = c(0, vapply(th, function(v) mean(p$prob[l] >= v), 1)))
  auc <- sum(diff(d$fpr) * (d$tpr[-1L] + d$tpr[-nrow(d)]) / 2)
  .plines(d$fpr, list(model = d$tpr),
          .pmain(main, sprintf("Linkage ROC curve (area %.3f)", auc)),
          "Share of unlinked pairs above the threshold",
          "Share of linked pairs above the threshold", col, pal, ...)
  graphics::abline(0, 1, lty = 3, col = "#636363")
  attr(d, "auc") <- auc
  invisible(d)
}

#' @rdname plotm
#' @export
plot.cconc <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cconc", c("lorenz", "stat"))
  g <- x$stat[["gini"]]
  d <- plotlrnz(x$lorenz, ...,
                main = .pmain(main, sprintf("Lorenz curve (Gini %.3f)", g)),
                col = col, pal = pal)
  attr(d, "gini") <- g
  invisible(d)
}

#' @rdname plotm
#' @export
plot.crossk <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "crossk", "k")
  k <- x$k
  .plines(k$r, list(estimated = k$k, `under independence` = k$theo),
          .pmain(main, sprintf("Cross K function, %s to %s", x$i, x$j)),
          "Distance", "K(r)", col, pal, where = "topleft", ...)
  invisible(k)
}

#' @rdname plotm
#' @export
plot.bandx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "bandx", "sens")
  s <- x$sens
  .plines(s$cap, list(value = s$value),
          .pmain(main, "The statistic by the cap on the open top band"),
          "Cap", "Statistic", col, pal, ...)
  invisible(s)
}

#' @rdname plotm
#' @export
plot.jdgiv <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "jdgiv", "sens")
  s <- x$sens
  .pgrp(s$dsh, s$cshare, paste("direct effect", s$deff),
        .pmain(main, "Complier share by the share of defiers"),
        "Share of defiers", "Complier share", col, pal, ...)
  invisible(s)
}

#' @rdname plotm
#' @export
plot.xcser <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "xcser", c("ccf", "cor"))
  c0 <- x$ccf
  b <- stats::qnorm(1 - (1 - x$level) / 2) / sqrt(x$cor$n[1L])
  cl <- .pcols(1L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(c0$lag), .prange(c0$r, b, -b),
          .pmain(main, "Cross-correlation at lags"), "Lag", "Correlation",
          ...)
  graphics::abline(h = 0, col = "#636363")
  graphics::abline(h = c(-b, b), lty = 2, col = "#636363")
  graphics::segments(c0$lag, 0, c0$lag, c0$r, col = cl, lwd = 3)
  graphics::points(c0$lag, c0$r, pch = 19, col = cl)
  invisible(c0)
}

#' @rdname plotm
#' @export
plot.hwkres <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "hwkres", c("u", "n"))
  n <- x$n
  q <- (seq_len(n) - 0.5) / n
  b <- 1.36 / sqrt(n)
  d <- data.frame(q = q, u = sort(x$u), lo = q - b, hi = q + b)
  cl <- .pcols(1L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(c(0, 1), c(0, 1),
          .pmain(main, "Time-rescaling residuals (95% Kolmogorov band)"),
          "uniform quantile", "rescaled time quantile", ...)
  graphics::abline(0, 1)
  graphics::abline(b, 1, lty = 2)
  graphics::abline(-b, 1, lty = 2)
  graphics::points(d$q, d$u, col = cl, pch = 1)
  invisible(d)
}

# ---- bars -------------------------------------------------------------------

#' @rdname plotm
#' @export
plot.arwt <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "arwt", "out")
  o <- x$out
  .pbars(o$value, o$to, .pmain(main, "Values moved to each target area"),
         "Value", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.crsev <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "crsev", "tab")
  r <- x$tab
  .pbars(r$idx, as.character(r$unit),
         .pmain(main, "Severity-weighted crime index"), "Index", col, pal,
         ...)
}

#' @rdname plotm
#' @export
plot.gmixlr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "gmixlr", "locus")
  r <- x$locus
  .pbars(r$log10lr, r$locus,
         .pmain(main, sprintf("log10 likelihood ratio by locus (total %.2f)",
                              x$log10lr)),
         "log10 likelihood ratio", col, pal, ref = 0, ...)
}

#' @rdname plotm
#' @export
plot.pedlr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "pedlr", "locus")
  r <- x$locus
  d <- .pbars(log10(r$lr), r$locus,
              .pmain(main, "log10 likelihood ratio by locus"),
              "log10 likelihood ratio", col, pal, ref = 0, ...)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.hypq <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "hypq", "wl")
  .pbars(x$wl, paste("unit", seq_along(x$wl)),
         .pmain(main, "Workload of each unit"), "Share of time busy", col,
         pal, ...)
}

#' @rdname plotm
#' @export
plot.lnsize <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "lnsize", "count")
  lab <- paste("member", seq_along(x$count))
  if (length(x$susp) == 1L) lab[x$susp] <- paste(lab[x$susp], "(suspect)")
  .pbars(x$count, lab,
         .pmain(main, "Choices of mock witnesses (line: an even spread)"),
         "Choices", col, pal, ref = mean(x$count), ...)
}

#' @rdname plotm
#' @export
plot.lostay <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "lostay", "bands")
  b <- x$bands
  lab <- ifelse(is.finite(b$hi), sprintf("%g to %g", b$lo, b$hi),
                sprintf("%g or more", b$lo))
  .pbars(b$share, lab, .pmain(main, "Lengths of stay"), "Share of stays",
         col, pal, ...)
}

#' @rdname plotm
#' @export
plot.pinpl <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "pinpl", "regs")
  r <- x$regs
  .pbars(r$units, as.character(r$region),
         .pmain(main, "Units in each region"), "Units", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.mmshr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "mmshr", c("group", "below", "at", "above"))
  d <- data.frame(group = x$group, below = x$below, at = x$at,
                  above = x$above, stringsAsFactors = FALSE)
  cl <- .pcols(3L, col, pal)
  op <- .ppar(mar = c(6, .plmar(d$group), 3, 1))
  on.exit(graphics::par(op), add = TRUE)
  m <- t(as.matrix(d[, -1L]))
  dn <- c(-1, 25, 12)
  graphics::barplot(m[, rev(seq_len(ncol(m))), drop = FALSE], horiz = TRUE,
                    names.arg = rev(d$group), col = cl, density = dn,
                    angle = c(0, 45, 135), border = "#252525",
                    xlab = "Share of sentences", main = "", ...)
  graphics::title(main = .pmain(main, "Sentences against the minimum"),
                  font.main = 1, cex.main = 1.05)
  graphics::legend("bottom", inset = c(0, -0.32), xpd = NA, horiz = TRUE,
                   legend = rownames(m), fill = cl, density = dn,
                   angle = c(0, 45, 135), bty = "n", cex = 0.85)
  invisible(d)
}

#' @rdname plotm
#' @export
plot.rptvic <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "rptvic", "tab")
  r <- x$tab
  .plines(r$k, list(observed = r$victims, `expected by chance` = r$expected),
          .pmain(main, "Targets by times victimised"), "Times victimised",
          "Targets", col, pal, ...)
  invisible(r[, c("k", "victims", "expected")])
}

# ---- matrices and surfaces ---------------------------------------------------

#' @rdname plotm
#' @export
plot.ctabx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "ctabx", c("tab", "res"))
  r <- x$res
  dimnames(r) <- dimnames(x$tab)
  .pmat(r, .pmain(main, "Residuals of the table (cells show counts)"),
        div = TRUE, lab = x$tab, pal = pal, ...)
}

#' @rdname plotm
#' @export
plot.knxnr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "knxnr", "tab")
  t <- x$tab
  rs <- paste(t$dlo, "to", t$dhi)
  cs <- paste(t$tlo, "to", t$thi)
  m <- tapply(t$ratio, list(factor(rs, unique(rs)), factor(cs, unique(cs))),
              identity)
  lab <- tapply(sprintf("%.2f%s", t$ratio, ifelse(t$p < 0.05, "*", "")),
                list(factor(rs, unique(rs)), factor(cs, unique(cs))),
                identity)
  .pmat(m, .pmain(main, "Knox ratios of close pairs (* p < 0.05)"),
        "Time apart", "Distance apart", div = TRUE, mid = 1, lab = lab,
        pal = pal, ...)
}

#' @rdname plotm
#' @export
plot.qdrtst <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "qdrtst", "counts")
  .pmat(x$counts, .pmain(main, "Counts by quadrat"), "Column", "Row",
        pal = pal, ...)
}

#' @rdname plotm
#' @export
plot.trnmat <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "trnmat", "prob")
  .pmat(x$prob, .pmain(main, "Transition probabilities, row to column"),
        "To", "From", pal = pal, ...)
}

#' @rdname plotm
#' @export
plot.zcomp <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "zcomp", "z")
  .pmat(x$z, .pmain(main, "Standardised counts by unit and category"),
        "Category", "Unit", div = TRUE, pal = pal, ...)
}

#' @rdname plotm
#' @export
plot.kdsurf <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "kdsurf", c("x", "y", "z"))
  .psurf(x$x, x$y, x$z, .pmain(main, "Kernel intensity surface"), pal = pal,
         ...)
}

#' @rdname plotm
#' @export
plot.prohs <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "prohs", c("x", "y", "z", "hot"))
  d <- .psurf(x$x, x$y, x$z,
              .pmain(main, "Prospective risk (hot cells outlined)"),
              pal = pal, ...)
  hx <- diff(x$x[1:2]) / 2
  hy <- diff(x$y[1:2]) / 2
  w <- which(x$hot, arr.ind = TRUE)
  graphics::rect(x$x[w[, 1L]] - hx, x$y[w[, 2L]] - hy, x$x[w[, 1L]] + hx,
                 x$y[w[, 2L]] + hy, border = "#252525", lwd = 1.5)
  d$hot <- x$hot
  invisible(d)
}

#' @rdname plotm
#' @export
plot.sdbsim <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "sdbsim", "b")
  .psurf(seq_len(nrow(x$b)), seq_len(ncol(x$b)), x$b,
         .pmain(main, "Attractiveness at the end of the run"), "column",
         "row", pal = pal, ...)
}

# ---- points -----------------------------------------------------------------

#' @rdname plotm
#' @export
plot.bufas <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "bufas", "pts")
  p <- x$pts
  g <- factor(p$grp, levels = unique(p$grp))
  cl <- .pcols(nlevels(g), col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(p$dsite), .prange(p$dcomp),
          .pmain(main, "Places by distance to sites and comparison sites"),
          "Distance to the nearest site",
          "Distance to the nearest comparison site", ...)
  graphics::abline(v = x$par, lty = c(2, 3), col = "#636363")
  graphics::points(p$dsite, p$dcomp, col = cl[g], pch = .ppch(g))
  .pkey("topright", levels(g), cl, pch = .ppch(seq_len(nlevels(g))))
  invisible(p)
}

#' @rdname plotm
#' @export
plot.cmpsc <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cmpsc", "segs")
  s <- x$segs
  cl <- .pcols(2L, col, pal)
  k <- s$cmps & !is.na(s$cmps)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(s$seg), .prange(s$ccp, x$pos),
          .pmain(main, sprintf("Segments congruent at the best shift: %d of %d",
                               x$score, x$nseg)),
          "Segment", "Position of the best match", ...)
  graphics::abline(h = x$pos, lty = 2, col = "#636363")
  graphics::points(s$seg, s$ccp, pch = ifelse(k, 19, 1), col = cl[k + 1L],
                   cex = 1.3)
  .pkey("topright", c("congruent", "not congruent"), cl[2:1], pch = c(19, 1))
  invisible(s)
}

#' @rdname plotm
#' @export
plot.ebshr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "ebshr", "est")
  e <- x$est
  cl <- .pcols(2L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(c(0.8, 2.2), .prange(e$y, e$shrunk),
          .pmain(main, "Counts and their empirical-Bayes shrinkage"), "",
          "Value", xaxt = "n", ...)
  graphics::axis(1, at = 1:2, labels = c("observed", "shrunk"),
                 tick = FALSE)
  graphics::segments(1, e$y, 2, e$shrunk, col = "#969696")
  graphics::points(rep(1, nrow(e)), e$y, pch = 1, col = cl[1L], cex = 1.2)
  graphics::points(rep(2, nrow(e)), e$shrunk, pch = 19, col = cl[2L],
                   cex = 1.2)
  invisible(e)
}

#' @rdname plotm
#' @export
plot.ibench <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "ibench", "tab")
  t <- x$tab
  f <- t$flag & !is.na(t$flag)
  cl <- .pcols(2L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  r <- .prange(t$bench, t$rate)
  .pfrm(r, r, .pmain(main, "Officer rate against benchmark"),
          "Benchmark", "Rate", ...)
  graphics::abline(0, 1, lty = 2, col = "#636363")
  graphics::points(t$bench, t$rate, pch = ifelse(f, 17, 1), col = cl[f + 1L],
                   cex = 1.2)
  if (any(f)) {
    graphics::text(t$bench[f], t$rate[f], t$off[f], pos = 3, cex = 0.75)
  }
  .pkey("topleft", c("flagged", "not flagged"), cl[2:1], pch = c(17, 1))
  invisible(t[, c("off", "bench", "rate", "flag")])
}

#' @rdname plotm
#' @export
plot.kscan <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "kscan", "cl")
  k <- x$cl
  cl <- .pcols(1L, col, pal)
  op <- .ppar()
  on.exit(graphics::par(op), add = TRUE)
  .pfrm(.prange(k$x - k$r, k$x + k$r), .prange(k$y - k$r, k$y + k$r),
          .pmain(main, "Space-time clusters (numbered by rank)"), "x", "y",
          asp = 1, ...)
  graphics::symbols(k$x, k$y, circles = k$r, inches = FALSE, add = TRUE,
                    fg = cl, lwd = 2)
  graphics::text(k$x, k$y, seq_len(nrow(k)), cex = 0.9)
  invisible(k)
}

# the event-time estimates of an event study or a staggered design
.pev <- function(x) {
  d <- if (inherits(x, "stgdid")) x$dynamic else x$est
  d[, c("evtm", "est", "lo", "hi")]
}
