#' Permutation importance of a prediction function
#'
#' How much worse a risk tool predicts when one of its inputs is made
#' useless: the column is shuffled across the rows, which keeps its
#' distribution but breaks its link with the outcome and with the other
#' columns, and the loss of the predictions against `y` is measured again
#' (Breiman 2001; Fisher, Rudin and Dominici 2019, who call the ratio the
#' model reliance). Each column is shuffled `nrep` times and the mean and
#' the standard deviation of the increases are reported.
#'
#' `f` is the prediction function: an R function that takes a numeric
#' matrix with the columns of `x` and returns one prediction per row, or
#' the names of learners of the package, fitted to `x` and `y` by
#' [suplrn()] with the settings in `lctl` (one learner is that learner, as
#' [rfor()] or least squares; several are their stack). The function is
#' called once, on the rows of `x` followed by every shuffled copy, so a
#' learner is fitted once. With learners the loss is measured on the rows
#' they were fitted to; for the loss on new rows, give `f` as a function
#' and `x` and `y` as the new rows.
#'
#' The shuffles are drawn in the compiled core from `seed`, the same on
#' every platform; `perm` returns them, so that every value can be checked.
#'
#' @inheritParams cemtch
#' @param x The rows at which to predict: a numeric matrix, its column
#'   names those `f` sees, no missing values.
#' @param f The prediction function: a function of a numeric matrix
#'   returning one finite number per row, or the names of learners of
#'   [suplrn()] (`"glm"`, `"lasso"`, `"ridge"`, `"enet"`, `"forest"`).
#' @param y The outcome at the rows of `x`: one finite number per row.
#' @param nrep The shuffles of each column: one whole number from 1.
#' @param loss `"se"`, the squared error (the Brier score for a 0/1
#'   outcome), or `"ae"`, the absolute error.
#' @param ratio `FALSE` for the increase of the loss, the loss after
#'   shuffling less the loss before; `TRUE` for their ratio.
#' @param seed The seed of the shuffles: one number, at least 0.
#' @param lctl The settings of the learners named in `f`: a list named by
#'   arguments of [suplrn()] (`fam`, `nfold`, `fold`, `disc`, `alpha`,
#'   `nlam`, `nfin`, `ntree`, `mtry`, `seed`); empty for a function.
#' @return An object of class `xpimp`, a list: `imp`, a data frame with
#'   one row per column of `x` (`term`, `imp` the mean increase or ratio,
#'   `sd` its standard deviation over the shuffles); `base`, the loss
#'   before shuffling; `rep`, the increase or ratio of every shuffle (one
#'   row per column, one column per shuffle); `perm`, the shuffles (an
#'   array of `nrow(x)` rows, `nrep` columns and one slice per column of
#'   `x`: `perm[i, r, j]` is the row whose value of column `j` row `i`
#'   takes in shuffle `r`); `loss` and `ratio`.
#' @references
#'   Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
#'   \doi{10.1023/A:1010933404324}
#'
#'   Fisher, A., Rudin, C. and Dominici, F. (2019). All models are wrong,
#'   but many are useful: learning a variable's importance by studying an
#'   entire class of prediction models simultaneously. *Journal of Machine
#'   Learning Research*, 20(177), 1-81.
#' @seealso [xppdp()], [xpale()], [xpcp()] and [xpshap()] for the shape of
#'   an input's effect; [fair()] and [calaud()] for the audit of a risk
#'   tool's scores.
#' @examples
#' set.seed(1)
#' x <- cbind(prior = rpois(200, 2), age = runif(200, 18, 60),
#'            noise = rnorm(200))
#' y <- 0.3 * x[, "prior"] - 0.02 * x[, "age"] + rnorm(200, sd = 0.5)
#' f <- function(m) 0.3 * m[, "prior"] - 0.02 * m[, "age"]
#' xpimp(x, f, y, nrep = 5)
#' xpimp(x, "glm", y, nrep = 5, lctl = list(nfold = 2))
#' @export
xpimp <- function(x, f, y, nrep = 10, loss = c("se", "ae"), ratio = FALSE,
                  seed = 1, lctl = list()) {
  fn <- "xpimp"
  a <- .xpck(x, f, y, lctl, fn)
  x <- a$x
  n <- nrow(x)
  p <- ncol(x)
  y <- .xpy(y, n, fn)
  nrep <- .cknum(nrep, "nrep", fn, min = 1, whole = TRUE)
  loss <- .ckarg(loss, c("se", "ae"), "loss", fn)
  if (!isTRUE(ratio) && !isFALSE(ratio)) {
    .ckstop(fn, "`ratio` must be TRUE or FALSE")
  }
  seed <- .cknum(seed, "seed", fn, min = 0)
  .xpbig(n * p * nrep, p, "nrep", fn)
  d <- .Call(C_clab_r_xpperm, x, nrep, seed)
  v <- .xpf(f, x, y, lctl, rbind(x, d$des), fn)
  k <- .Call(C_clab_r_xpimp, y, v[seq_len(n)], v[-seq_len(n)], p, nrep,
             match(loss, c("se", "ae")) - 1L, ratio)
  rep <- k$rep
  dimnames(rep) <- list(a$nm, NULL)
  st <- k$st
  structure(
    list(imp = data.frame(term = a$nm, imp = st[1L + seq_len(p)],
                          sd = st[1L + p + seq_len(p)]),
         base = st[1L], rep = rep, perm = d$perm, loss = loss,
         ratio = ratio),
    class = "xpimp"
  )
}

#' @rdname xpimp
#' @export
print.xpimp <- function(x, ...) {
  cat("xpimp: permutation importance, ", ncol(x[["rep"]]),
      " shuffles per column; loss ", x[["loss"]],
      if (x[["ratio"]]) " (ratio)" else " (increase)", ", ",
      format(x[["base"]], digits = 4), " before shuffling\n", sep = "")
  i <- x[["imp"]]
  print(i[order(-i[["imp"]]), , drop = FALSE], row.names = FALSE)
  invisible(x)
}

#' Partial dependence of a prediction function on one input
#'
#' The average prediction when one input is set to a value for every row
#' and the others keep their values (Friedman 2001, sec. 8.2): at each
#' value `g` of a grid, the mean of the predictions with column `v` set to
#' `g` in every row. It is
#' the effect of the input averaged over the data, as a risk tool's user
#' would read it; where the input is correlated with the others it averages
#' over combinations the data never hold, which [xpale()] avoids.
#'
#' The grid is `grid` values of the column, evenly spaced from its minimum
#' to its maximum (`gtype = "even"`) or at its quantiles
#' (`gtype = "quant"`, `quantile()`'s default type 7, repeated values
#' dropped), or the values given. The rows are written in the compiled
#' core; `f` is called once on all of them.
#'
#' @inheritParams xpimp
#' @param v The input: the name or the number of one column of `x`.
#' @param y The outcome at the rows of `x`, needed only to fit the learners
#'   named in `f`.
#' @param grid The number of grid points, one whole number from 2, or the
#'   values of the grid (at least two).
#' @param gtype `"even"` for points evenly spaced over the column's range,
#'   `"quant"` for its quantiles.
#' @return An object of class `xppdp`, a list: `pd`, a data frame with the
#'   grid (`grid`) and the partial dependence there (`pd`); `term`, the
#'   input's name; and `n`, the rows averaged.
#' @references
#'   Friedman, J. H. (2001). Greedy function approximation: a gradient
#'   boosting machine. *Annals of Statistics*, 29(5), 1189-1232.
#'   \doi{10.1214/aos/1013203451}
#' @seealso [xpcp()] for the same curve for single rows, [xpale()].
#' @examples
#' set.seed(2)
#' x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
#' f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
#' xppdp(x, f, "age", grid = 5)
#' @export
xppdp <- function(x, f, v, y = NULL, grid = 20, gtype = c("even", "quant"),
                  lctl = list()) {
  fn <- "xppdp"
  a <- .xpck(x, f, y, lctl, fn)
  x <- a$x
  j <- .xpv(v, a$nm, fn)
  g <- .xpg(x, j, grid, gtype, fn)
  .xpbig(nrow(x) * length(g), ncol(x), "grid", fn)
  d <- .Call(C_clab_r_xpcpd, x, j, g)
  colnames(d) <- colnames(x)
  pd <- .Call(C_clab_r_xppd, .xpf(f, x, y, lctl, d, fn), nrow(x), length(g))
  structure(list(pd = data.frame(grid = g, pd = pd), term = a$nm[j],
                 n = nrow(x)),
            class = "xppdp")
}

#' @rdname xppdp
#' @export
print.xppdp <- function(x, ...) {
  cat("xppdp: partial dependence on ", x[["term"]], ", averaged over ",
      x[["n"]], " rows\n", sep = "")
  print(x[["pd"]], row.names = FALSE)
  invisible(x)
}

#' Ceteris paribus profiles of single rows
#'
#' The prediction for one person as one input moves over a grid while the
#' person's other inputs stay as they are: what the risk tool would have
#' said had only that input been different. One profile per row of `newx`;
#' over all the rows of `x` the profiles are the individual conditional
#' expectation curves of Goldstein et al. (2015), and their mean is the
#' partial dependence of [xppdp()]. A steep profile on a protected or
#' proxy input is a decision that turns on that input alone.
#'
#' The grid is taken from the column of `x`, as in [xppdp()]. The rows are
#' written in the compiled core; `f` is called once on all of them.
#'
#' @inheritParams xppdp
#' @param newx The rows to profile: a numeric matrix with the columns of
#'   `x` (by default its first row).
#' @return An object of class `xpcp`, a list: `prof`, the predictions (one
#'   row per row of `newx`, one column per grid value); `grid`; `pred`, the
#'   prediction at each row of `newx` as it is; and `term`.
#' @references
#'   Goldstein, A., Kapelner, A., Bleich, J. and Pitkin, E. (2015). Peeking
#'   inside the black box: visualizing statistical learning with plots of
#'   individual conditional expectation. *Journal of Computational and
#'   Graphical Statistics*, 24(1), 44-65.
#'   \doi{10.1080/10618600.2014.907095}
#'
#'   Biecek, P. and Burzykowski, T. (2021). *Explanatory Model Analysis*.
#'   Chapman and Hall/CRC, New York, ch. 10.
#' @examples
#' set.seed(3)
#' x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
#' f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
#' xpcp(x, f, "prior", newx = x[1:2, ], grid = c(0, 2, 4, 6))
#' @export
xpcp <- function(x, f, v, newx = x[1L, , drop = FALSE], y = NULL, grid = 20,
                 gtype = c("even", "quant"), lctl = list()) {
  fn <- "xpcp"
  a <- .xpck(x, f, y, lctl, fn)
  x <- a$x
  j <- .xpv(v, a$nm, fn)
  newx <- .xpnew(newx, x, fn)
  g <- .xpg(x, j, grid, gtype, fn)
  ns <- nrow(newx)
  .xpbig(ns * (length(g) + 1), ncol(x), "grid", fn)
  d <- .Call(C_clab_r_xpcpd, newx, j, g)
  colnames(d) <- colnames(x)
  k <- .xpf(f, x, y, lctl, rbind(d, newx), fn)
  prof <- matrix(k[seq_len(ns * length(g))], ns, length(g))
  structure(list(prof = prof, grid = g, pred = k[-seq_len(ns * length(g))],
                 term = a$nm[j]),
            class = "xpcp")
}

#' @rdname xpcp
#' @export
print.xpcp <- function(x, ...) {
  r <- apply(x[["prof"]], 1L, range)
  cat("xpcp: ceteris paribus profiles of ", nrow(x[["prof"]]),
      " row(s) over ", length(x[["grid"]]), " values of ", x[["term"]],
      "\n", sep = "")
  print(data.frame(pred = x[["pred"]], lo = r[1L, ], hi = r[2L, ]),
        row.names = FALSE)
  invisible(x)
}

#' Accumulated local effects of one input
#'
#' The effect of one input on a prediction, measured only where the data
#' are (Apley and Zhu 2020): the range of the column is cut at its
#' quantiles into intervals, each row's prediction is taken at the upper
#' and at the lower edge of its own interval with its other inputs as they
#' are, the differences are averaged within each interval and accumulated
#' from the lowest edge. Unlike partial dependence ([xppdp()]) no row is
#' moved far from where it lies, so an input correlated with others (prior
#' arrests and age, say) is not credited with their effects.
#'
#' The edges are the column's minimum and its quantiles at
#' `1/nbin, 2/nbin, ..., 1` (`quantile()` type 1), repeated values dropped,
#' so there are at most `nbin` intervals; each interval holds the rows
#' above its lower edge up to its upper one, the first its lower edge too.
#' The effects are centred so that their mean over the rows, the effect
#' taken as linear within each interval, is 0: the centring of the
#' authors' 'ALEPlot'. The rows are written and the effects accumulated in
#' the compiled core; `f` is called once on all the rows.
#'
#' @inheritParams xppdp
#' @param nbin The intervals asked for: one whole number from 1 to
#'   67108864.
#' @return An object of class `xpale`, a list: `ale`, a data frame with
#'   the edges (`z`) and the centred accumulated effect there (`ale`);
#'   `bins`, a data frame with one row per interval (`lo`, `hi` its edges,
#'   `n` its rows, `eff` the mean difference of the predictions across it);
#'   and `term`.
#' @references
#'   Apley, D. W. and Zhu, J. (2020). Visualizing the effects of predictor
#'   variables in black box supervised learning models. *Journal of the
#'   Royal Statistical Society Series B*, 82(4), 1059-1086.
#'   \doi{10.1111/rssb.12377}
#' @examples
#' set.seed(4)
#' age <- runif(200, 18, 60)
#' prior <- rpois(200, (60 - age) / 10)
#' x <- cbind(prior, age)
#' f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
#' xpale(x, f, "age", nbin = 5)
#' @export
xpale <- function(x, f, v, y = NULL, nbin = 10, lctl = list()) {
  fn <- "xpale"
  a <- .xpck(x, f, y, lctl, fn)
  x <- a$x
  j <- .xpv(v, a$nm, fn)
  nbin <- .cknum(nbin, "nbin", fn, min = 1, whole = TRUE)
  if (nbin > 67108864) .ckstop(fn, "`nbin` must be at most 67108864")
  .xpbig(2 * nrow(x), ncol(x), "x", fn)
  d <- .Call(C_clab_r_xpaled, x, j, nbin)
  colnames(d$des) <- colnames(x)
  z <- d$z
  k <- .Call(C_clab_r_xpale, .xpf(f, x, y, lctl, d$des, fn), d$bin,
             length(z) - 1L)
  structure(
    list(ale = data.frame(z = z, ale = k$ale),
         bins = data.frame(lo = z[-length(z)], hi = z[-1L], n = k$cnt,
                           eff = k$eff),
         term = a$nm[j]),
    class = "xpale"
  )
}

#' @rdname xpale
#' @export
print.xpale <- function(x, ...) {
  cat("xpale: accumulated local effects of ", x[["term"]], " over ",
      nrow(x[["bins"]]), " interval(s)\n", sep = "")
  print(x[["ale"]], row.names = FALSE)
  invisible(x)
}

#' Sampling Shapley values of single predictions
#'
#' Shares out one person's prediction among the inputs: the Shapley value
#' of input `j` is its contribution averaged over every order in which the
#' inputs could be revealed, the inputs not yet revealed taking values from
#' the rows of `x` (the reference data). Strumbelj and Kononenko (2014)
#' estimate it by sampling: each sample draws an order of the inputs and a
#' row `w` of `x`, and takes the prediction at the person's values on `j`
#' and the inputs before it (`w`'s on the rest) less the prediction with
#' `j` too taken from `w`. The value is the mean over `nsamp` samples, its
#' standard error their standard deviation over `sqrt(nsamp)`. Over all the
#' inputs the values add up, in expectation, to the person's prediction
#' less the mean prediction over `x`.
#'
#' The orders and rows are drawn in the compiled core from `seed`, the same
#' on every platform; `f` is called once on all the rows.
#'
#' @inheritParams xpcp
#' @param x The reference data: a numeric matrix, its column names those
#'   `f` sees, no missing values (the rows the learners in `f` are fitted
#'   to).
#' @param newx The rows to explain: a numeric matrix with the columns of
#'   `x` (by default its first row).
#' @param nsamp The samples of each value: one whole number from 1.
#' @param seed The seed of the draws: one number, at least 0.
#' @return An object of class `xpshap`, a list: `phi` and `se`, the values
#'   and their standard errors (one row per row of `newx`, one column per
#'   input); `pred`, the prediction at each row of `newx`; `base`, the mean
#'   prediction over the rows of `x`; and `nsamp`.
#' @references
#'   Strumbelj, E. and Kononenko, I. (2014). Explaining prediction models
#'   and individual predictions with feature contributions. *Knowledge and
#'   Information Systems*, 41(3), 647-665. \doi{10.1007/s10115-013-0679-x}
#'
#'   Shapley, L. S. (1953). A value for n-person games. In Kuhn, H. W. and
#'   Tucker, A. W. (eds), *Contributions to the Theory of Games II*,
#'   307-317. Princeton University Press.
#' @examples
#' set.seed(5)
#' x <- cbind(prior = rpois(100, 2), age = runif(100, 18, 60))
#' f <- function(m) plogis(-1 + 0.5 * m[, "prior"] - 0.03 * m[, "age"])
#' xpshap(x, f, newx = x[1:2, ], nsamp = 100)
#' @export
xpshap <- function(x, f, newx = x[1L, , drop = FALSE], y = NULL,
                   nsamp = 200, seed = 1, lctl = list()) {
  fn <- "xpshap"
  a <- .xpck(x, f, y, lctl, fn)
  x <- a$x
  newx <- .xpnew(newx, x, fn)
  nsamp <- .cknum(nsamp, "nsamp", fn, min = 1, whole = TRUE)
  seed <- .cknum(seed, "seed", fn, min = 0)
  n <- nrow(x)
  ns <- nrow(newx)
  p <- ncol(x)
  .xpbig(2 * ns * p * nsamp + ns + n, p, "nsamp", fn)
  d <- .Call(C_clab_r_xpshd, x, newx, nsamp, seed)
  colnames(d) <- colnames(x)
  k <- .xpf(f, x, y, lctl, rbind(d, newx, x), fn)
  q <- nrow(d)
  s <- .Call(C_clab_r_xpsh, k[seq_len(q)], ns, p, nsamp)
  dn <- list(NULL, a$nm)
  b <- .Call(C_clab_r_xppd, k[q + ns + seq_len(n)], n, 1L)
  structure(
    list(phi = `dimnames<-`(s$phi, dn), se = `dimnames<-`(s$se, dn),
         pred = k[q + seq_len(ns)], base = b,
         nsamp = nsamp),
    class = "xpshap"
  )
}

#' @rdname xpshap
#' @export
print.xpshap <- function(x, ...) {
  cat("xpshap: sampling Shapley values, ", x[["nsamp"]],
      " samples each; mean prediction ", format(x[["base"]], digits = 4),
      "\n", sep = "")
  print(data.frame(pred = x[["pred"]], x[["phi"]], check.names = FALSE),
        row.names = FALSE)
  invisible(x)
}

# The settings of suplrn() that `lctl` may hold.
.xpset <- c("fam", "nfold", "fold", "disc", "alpha", "nlam", "nfin", "ntree",
            "mtry", "seed")

# The checks the explanation functions share: x a numeric matrix, finite,
# its column names (x1, x2, ... where it has none); f a function or the
# names of learners; lctl their settings; y given for learners.
.xpck <- function(x, f, y, lctl, fn) {
  if (!is.matrix(x) || !is.numeric(x) || !nrow(x) || !ncol(x)) {
    .ckstop(fn, paste("`x` must be a numeric matrix with at least one row",
                      "and one column"))
  }
  if (!all(is.finite(x))) {
    .ckstop(fn, "`x` must be finite numbers with no missing values")
  }
  storage.mode(x) <- "double"
  nm <- colnames(x)
  if (is.null(nm)) nm <- rep("", ncol(x))
  e <- is.na(nm) | !nzchar(nm)
  nm[e] <- paste0("x", which(e))
  lrn <- is.character(f) && length(f) > 0L && !anyNA(f)
  if (!is.function(f) && !lrn) {
    .ckstop(fn, "`f` must be a function or the names of learners of suplrn()")
  }
  k <- names(lctl)
  if (!is.list(lctl) || is.data.frame(lctl) ||
        (length(lctl) && (is.null(k) || !all(k %in% .xpset) ||
                            anyDuplicated(k)))) {
    .ckstop(fn, "`lctl` must be a list named by settings of suplrn(): %s",
            toString(.xpset))
  }
  if (!lrn && length(lctl)) {
    .ckstop(fn, "`lctl` is for learners named in `f`, not for a function")
  }
  if (lrn && is.null(y)) {
    .ckstop(fn, "`y` must be given to fit the learners named in `f`")
  }
  list(x = x, nm = nm)
}

# The outcome: one finite number per row of x.
.xpy <- function(y, n, fn) {
  y <- as.numeric(.cknum(y, "y", fn, scalar = FALSE))
  if (length(y) != n) {
    .ckstop(fn, "`y` must have one value per row of `x`, %d", n)
  }
  y
}

# A design of rows times p values within 2^26.
.xpbig <- function(rows, p, what, fn) {
  if (rows * p > 67108864) {
    .ckstop(fn, paste("`%s` asks for %.0f rows of %d columns to predict,",
                      "more than 67108864 values"), what, rows, p)
  }
}

# The predictions of f at the rows of m: f(m) for a function, else the
# learners named by f fitted to x and y by suplrn() with the settings lctl.
.xpf <- function(f, x, y, lctl, m, fn) {
  if (is.function(f)) {
    v <- f(m)
    if (!is.numeric(v) || length(v) != nrow(m) || !all(is.finite(v))) {
      .ckstop(fn, "`f` must return a finite number for each of the %d rows",
              nrow(m))
    }
    return(as.numeric(v))
  }
  a <- c(list(x = x, y = y, lrns = f, newx = m), lctl)
  .clas(fn, do.call(suplrn, a))$pred
}

# The column v of x, by name or number: its number.
.xpv <- function(v, nm, fn) {
  j <- if (is.character(v) && length(v) == 1L && !is.na(v)) {
    match(v, nm)
  } else if (is.numeric(v) && length(v) == 1L && is.finite(v) &&
               v == round(v) && v >= 1 && v <= length(nm)) {
    as.integer(v)
  } else {
    NA_integer_
  }
  if (is.na(j)) .ckstop(fn, "`v` must name or number one column of `x`")
  j
}

# The grid of column j of x: grid points spaced by gtype, or the values
# given.
.xpg <- function(x, j, grid, gtype, fn) {
  gtype <- .ckarg(gtype, c("even", "quant"), "gtype", fn)
  ok <- is.numeric(grid) && length(grid) && all(is.finite(grid))
  if (ok && length(grid) == 1L) {
    ok <- grid >= 2 && grid <= 67108864 && grid == round(grid)
  }
  if (!ok) {
    .ckstop(fn, paste("`grid` must be a whole number of points from 2, or",
                      "at least two values"))
  }
  if (length(grid) > 1L) return(as.numeric(grid))
  .Call(C_clab_r_xpgrid, x[, j], grid, match(gtype, c("even", "quant")) - 1L)
}

# The instances: a numeric matrix with the columns of x, finite.
.xpnew <- function(newx, x, fn) {
  if (!is.matrix(newx) || !is.numeric(newx) || !nrow(newx) ||
        ncol(newx) != ncol(x) || !all(is.finite(newx))) {
    .ckstop(fn, paste("`newx` must be a numeric matrix with the %d columns",
                      "of `x`, no missing values"), ncol(x))
  }
  storage.mode(newx) <- "double"
  colnames(newx) <- colnames(x)
  newx
}
