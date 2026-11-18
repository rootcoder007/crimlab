# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# Plot methods for the results on Indigenous health and wellbeing. Each one
# is a dot-and-interval chart of what the result's print method tabulates:
# a group's level beside its interval, never a ranking drawn without it.

# a group's label, starred where the source flags the estimate (E)
.pstar <- function(lab, flag) {
  lab <- as.character(lab)
  f <- !is.na(flag) & as.logical(flag)
  lab[f] <- paste(lab[f], "*")
  lab
}

#' @rdname plotm
#' @export
plot.indgap <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "indgap", c("group", "ratio"))
  pooled <- !is.na(x$poogrp) & x$poogrp
  g <- if (any(pooled)) {
    ifelse(pooled, "pooled total", "distinct group")
  }
  .pdots(.pstar(x$group, x$usewithcau), x$ratio, x$ratiolo, x$ratioup,
         .pmain(main, paste("Ratio to", attr(x, "reference"))),
         "Ratio (log scale); 1 is no gap", col, pal, ref = 1, grp = g,
         log = "x", ...)
}

#' @rdname plotm
#' @export
plot.dualcont <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "dualcont", c("group", "wellness", "distress"))
  lab <- attr(x, "labels")
  k <- if (length(lab) == 2L) {
    paste(c("wellness:", "distress:"), lab)
  } else {
    c("wellness", "distress")
  }
  .pdots(rep(.pstar(x$group, x$usewithcau), 2L), c(x$wellness, x$distress),
         c(x$wello, x$distrlo), c(x$welup, x$distrup),
         .pmain(main, "Wellness beside distress"), "Estimate", col, pal,
         grp = rep(k, each = nrow(x)), ...)
}

#' @rdname plotm
#' @export
plot.commidx <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "commidx", c("group", "index"))
  .pdots(.pstar(x$group, x$usewithcau), x$index, x$lower, x$upper,
         .pmain(main, "Community-weighted index"),
         "Index, under the weights given", col, pal, ...)
}

#' @rdname plotm
#' @export
plot.stdprev <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "stdprev", c("group", "estimate"))
  .pdots(.pstar(x$group, x$usewithcau), x$estimate, x$lower, x$upper,
         .pmain(main, "Age-standardised prevalence"), "Prevalence", col,
         pal, ...)
}

#' @rdname plotm
#' @export
plot.ordcmp <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "ordcmp", c("group", "probbet"))
  .pdots(.pstar(x$group, x$usewithcau), x$probbet, x$lower, x$upper,
         .pmain(main, paste("Chance of a higher rating than",
                            x$reference[1L])),
         "Probability; 0.5 is no difference", col, pal, ref = 0.5, ...)
}

#' @rdname plotm
#' @export
plot.welprof <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "welprof", c("group", "indicator", "estimate"))
  x <- x[!is.na(x$estimate), , drop = FALSE]
  .pdots(paste0(x$indicator, " (", x$side, ")"), x$estimate, x$lower,
         x$upper, .pmain(main, "Wellness profile"), "Estimate", col, pal,
         grp = .pstar(x$group, x$usewithcau), ...)
}

#' @rdname plotm
#' @export
plot.holbal <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "holbal", c("group", "level", "balance"))
  op <- .ppanels(2L)
  on.exit(graphics::par(op), add = TRUE)
  g <- .pstar(x$group, x$usewithcau)
  out <- list(level = .pdots(g, x$level, x$lvllo, x$lvlup,
                             .pmain(main, "Level across the domains"),
                             "Mean level", col, pal, ...))
  out$balance <- .pdots(g, x$balance, x$ballo, x$balup,
                        "Balance of the domains",
                        "Balance; 1 is even", col, pal,
                        ref = 1, ...)
  invisible(out)
}

#' @rdname plotm
#' @export
plot.healchg <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "healchg", c("group", "indicator", "change"))
  p <- attr(x, "periods")
  .pdots(paste0(x$indicator, " (", x$side, ")"), x$change, x$lower,
         x$upper,
         .pmain(main, paste("Change,", paste(p, collapse = " to "))),
         "Change in the estimate's units; 0 is none", col, pal, ref = 0,
         grp = .pstar(x$group, x$usewithcau), ...)
}

#' @rdname plotm
#' @export
plot.trchlthgap <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "trchlthgap", c("group", "disratio"))
  .pdots(.pstar(x$group, x$usewithcau), x$disratio, x$lower, x$upper,
         .pmain(main, "Health gaps (Call to Action 19)"),
         "Disparity ratio (log scale); above 1 is worse", col, pal,
         ref = 1, grp = paste(x$indicator, x$period), log = "x", ...)
}

# the four components of a Community Well-Being score and the score, one
# row per community or group
.pcwb <- function(lab, x, main, col, pal, ...) {
  k <- c(income = "income", education = "education", housing = "housing",
         labour = "labour force", cwb = "score")
  v <- unlist(x[names(k)], use.names = FALSE)
  .pdots(rep(lab, length(k)), v, NA, NA, main,
         "0 to 100; the score is the mean of the four", col, pal,
         grp = rep(k, each = length(lab)), ...)
}

#' @rdname plotm
#' @export
plot.cnscwb <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cnscwb", c("geo", "group", "income"))
  .pcwb(paste(x$geo, x$group, sep = ": "), x,
        .pmain(main, "Community Well-Being (census)"), col, pal, ...)
}

#' @rdname plotm
#' @export
plot.cwbscr <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "cwbscr", c("income", "cwb"))
  lab <- rownames(x)
  if (is.null(lab)) lab <- as.character(seq_len(nrow(x)))
  .pcwb(paste("community", lab), x,
        .pmain(main, "Community Well-Being scores"), col, pal, ...)
}
