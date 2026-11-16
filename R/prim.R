# Summary statistics computed by the compiled core (src/prim.cpp). The R
# functions check their arguments; .pcall() runs the kernel.

#' Mean, variance and correlation (C backend)
#'
#' `cmn()`, `cvari()` and `ccorr()` compute the mean, the sample variance
#' and the Pearson correlation in compiled code, the same kernels that
#' crimlab's other front ends (the command line, Python and the bindings)
#' call, so a result computed here is the result they give. For everyday
#' work `mean()`, `var()` and `cor()` give the same numbers.
#'
#' With \eqn{n} values,
#' \deqn{\bar x = \frac{1}{n}\sum_i x_i, \qquad
#'   s^2 = \frac{1}{n - 1}\sum_i (x_i - \bar x)^2, \qquad
#'   r = \frac{\sum_i (x_i - \bar x)(y_i - \bar y)}
#'   {\sqrt{\sum_i (x_i - \bar x)^2 \sum_i (y_i - \bar y)^2}}.}
#' The variance uses the \eqn{n - 1} (sample) denominator, as
#' [stats::var()] does. The mean is base R's (an extended-precision sum
#' and one corrective pass), and a sum of finite values that would
#' overflow a double (two values of `1e308`) still gives the finite mean.
#' The correlation is computed from centred sums and clamped to
#' \eqn{[-1, 1]}.
#'
#' There is no `na.rm`: an `NA` or `NaN` gives `NA` or `NaN`; remove
#' missing values first, for example with [stats::na.omit()]. Edge cases
#' return `NaN` rather than an error: an empty vector for the mean, fewer
#' than two values for the variance, a constant vector for the
#' correlation. Text, a factor, a list or a data frame is refused; logical
#' input is read as 0 and 1. `ccorr()` refuses `x` and `y` of different
#' lengths.
#'
#' @param x,y Numeric (or logical) vectors; `x` and `y` of the same length
#'   in `ccorr()`.
#' @return A number.
#' @seealso [cmom()], [cspr()], [ccovm()].
#' @examples
#' cmn(1:10)
#' cmn(c(1, 2, NA))
#' cmn(c(1e308, 1e308))
#' cvari(c(2, 4, 4, 4, 5, 5, 7, 9))
#' all.equal(cvari(1:10), stats::var(1:10))
#' ccorr(1:10, (1:10)^2)
#' ccorr(1:10, 10:1)
#' try(cmn(c("1", "2")))
#' try(ccorr(1:3, 1:4))
#' @name cstat
#' @export
cmn <- function(x) .pcall(0L, .pnum(x, "x", "cmn"))

#' @rdname cstat
#' @export
cvari <- function(x) .pcall(1L, .pnum(x, "x", "cvari"), a = 1)

#' @rdname cstat
#' @export
ccorr <- function(x, y) {
  x <- .pnum(x, "x", "ccorr")
  y <- .pnum(y, "y", "ccorr")
  .psame(x, y, "x", "y", "ccorr")
  .pcall(2L, x, y)
}

#' Standard deviation and Euclidean distance (C backend)
#'
#' `csdev()` returns the standard deviation of a numeric vector and
#' `cdistm()` the straight-line (Euclidean) distance between two vectors
#' of the same length.
#'
#' `csdev()` computes \deqn{s = \sqrt{\frac{1}{n - \mathrm{ddof}}
#' \sum_{i=1}^n (x_i - \bar x)^2},} so `ddof = 1` is the sample standard
#' deviation of [stats::sd()] and `ddof = 0` the population figure.
#' `cdistm()` computes \eqn{\sqrt{\sum_i (a_i - b_i)^2}}.
#'
#' Missing values propagate: one `NA` makes the result `NA`. A vector with
#' no more than `ddof` values has no defined spread and gives `NaN`. Two
#' empty vectors are at distance 0. Text, a factor, a list or a data frame
#' is refused, and so are `a` and `b` of different lengths; logical input
#' counts as 0 and 1.
#'
#' @param x,a,b Numeric (or logical) vectors; `a` and `b` of the same
#'   length.
#' @param ddof The delta degrees of freedom, one whole number from 0: the
#'   divisor is `n - ddof`.
#' @return A number: the standard deviation in the units of `x`, or the
#'   distance in the units of `a` and `b`.
#' @seealso [cmom()], [cmad()].
#' @examples
#' csdev(c(2, 4, 4, 4, 5, 5, 7, 9))
#' all.equal(csdev(1:10), stats::sd(1:10))
#' csdev(1:10, ddof = 0)
#' csdev(c(3, NA, 5))
#' cdistm(c(0, 0), c(3, 4))
#' try(cdistm(1:3, 1:2))
#' @name cspread
#' @export
csdev <- function(x, ddof = 1L) {
  x <- .pnum(x, "x", "csdev")
  d <- .pone(ddof, "ddof", "csdev", 0, .Machine$integer.max, TRUE)
  sqrt(.pcall(1L, x, a = d))
}

#' @rdname cspread
#' @export
cdistm <- function(a, b) {
  a <- .pnum(a, "a", "cdistm")
  b <- .pnum(b, "b", "cdistm")
  .psame(a, b, "a", "b", "cdistm")
  .pcall(3L, a, b)
}

#' Normal density and log-density (C backend)
#'
#' `cnormpdf()` evaluates the density of the normal distribution with a
#' given mean and standard deviation at each value of `x`; `cnormlpdf()`
#' its logarithm, computed directly from the formula rather than as
#' `log(dnorm(x))`. Summing log-densities is how a normal model is scored,
#' and computing them directly keeps the sum finite far into the tails,
#' where the density itself underflows to zero.
#'
#' \deqn{f(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x -
#' \mu)^2}{2\sigma^2}\right), \qquad \log f(x) = -\log \sigma -
#' \tfrac12 \log(2\pi) - \frac{(x - \mu)^2}{2\sigma^2}.}
#' Both agree with [stats::dnorm()] to rounding. They are vectorised over
#' `x` only: `mean` must be one finite number and `sd` one positive
#' finite number. `NA` in `x` gives `NA`; text, a factor or a list is
#' refused.
#'
#' @param x Numeric vector of quantiles.
#' @param mean The mean \eqn{\mu}, one finite number.
#' @param sd The standard deviation \eqn{\sigma}, one positive finite
#'   number.
#' @return A numeric vector the length of `x`.
#' @seealso [stats::dnorm()].
#' @examples
#' cnormpdf(c(-1, 0, 1))
#' cnormpdf(c(0, 5, 10), mean = 5, sd = 2)
#' all.equal(cnormpdf(-2:2), stats::dnorm(-2:2))
#' stats::dnorm(50)
#' cnormlpdf(50)
#' y <- c(4.8, 5.1, 5.6, 4.9, 5.3)
#' sum(cnormlpdf(y, mean = mean(y), sd = stats::sd(y)))
#' try(cnormpdf(0, sd = 0))
#' try(cnormlpdf(1, mean = c(0, 1)))
#' @name cnorm
#' @export
cnormpdf <- function(x, mean = 0, sd = 1) .pnorm(x, mean, sd, 0, "cnormpdf")

#' @rdname cnorm
#' @export
cnormlpdf <- function(x, mean = 0, sd = 1) .pnorm(x, mean, sd, 1, "cnormlpdf")

#' Mean, variance, skewness and kurtosis (C backend)
#'
#' The four standard descriptive moments of a numeric vector in one call:
#' where it is centred, how spread out it is, whether it leans to one side
#' and how heavy its tails are.
#'
#' The mean is base R's (an extended-precision sum and a corrective pass);
#' the deviations from it are divided by their largest magnitude before
#' any power is taken, and the scale is restored at the end, which avoids
#' both the cancellation of the textbook formula and the overflow of
#' fourth powers beyond about \eqn{10^{77}}. The variance uses the `n - 1`
#' denominator of [stats::var()]. The shape statistics are the sample
#' moment ratios \deqn{g_1 = m_3 / m_2^{3/2}, \qquad g_2 = m_4 / m_2^2 - 3,}
#' with \eqn{m_k = n^{-1} \sum_i (x_i - \bar x)^k}, so kurtosis is EXCESS
#' kurtosis (a normal sample sits near zero). These are the uncorrected
#' ratios (type 1 of Joanes and Gill, 1998).
#'
#' Edge cases are `NaN`, never an error: an empty vector gives four `NaN`;
#' one value gives its mean and `NaN` for the rest; skewness needs three
#' values and kurtosis four; any `NA` or `NaN` makes every entry `NaN`. A
#' constant vector has variance 0 and `NaN` shape statistics.
#'
#' @param x Numeric (or logical) vector.
#' @return A named numeric vector: `mean`, `variance`, `skewness` and
#'   `kurtosis`.
#' @references Joanes DN, Gill CA (1998). Comparing measures of sample
#'   skewness and kurtosis. *Journal of the Royal Statistical Society:
#'   Series D (The Statistician)* 47(1), 183--189.
#'   \doi{10.1111/1467-9884.00122}
#' @seealso [csdev()], [cqnt()], [cmad()].
#' @examples
#' cmom(c(2, 4, 4, 4, 5, 5, 7, 9))
#' m <- cmom(1:10)
#' all.equal(m[["variance"]], stats::var(1:10))
#' cmom(c(-2, -1, 0, 1, 2))[["skewness"]]
#' cmom(c(1, 2))
#' @export
cmom <- function(x) {
  out <- .pcall(5L, .pnum(x, "x", "cmom"))
  names(out) <- c("mean", "variance", "skewness", "kurtosis")
  out
}

#' Quantiles, median and robust spread (C backend)
#'
#' `cqnt()` is the type-7 quantile, R's default, so it agrees with
#' `stats::quantile(x, probs)`. `cmed()` is the 50% point. `cmad()` is the
#' median absolute deviation scaled by `const` so that it estimates the
#' standard deviation of a normal sample. `ciqr()` is the interquartile
#' range and `ctkyfnc()` the outlier fences drawn `k` IQRs beyond the
#' quartiles. A single corrupted row can move a mean or a variance
#' arbitrarily far, but moves a median or a MAD hardly at all, which is
#' what is wanted when screening a column of a released table for values
#' to inspect by hand.
#'
#' The type-7 quantile at probability \eqn{p} interpolates between the
#' order statistics \eqn{x_{(j)}} and \eqn{x_{(j+1)}} with \eqn{h = (n -
#' 1)p + 1}, \eqn{j = \lfloor h \rfloor}: \eqn{Q(p) = x_{(j)} + (h - j)
#' (x_{(j+1)} - x_{(j)})} (Hyndman and Fan, 1996). The MAD is \eqn{c \cdot
#' \mathrm{median}_i |x_i - \mathrm{median}(x)|}; the IQR \eqn{Q(0.75) -
#' Q(0.25)}; the fences \eqn{Q(0.25) - k \cdot \mathrm{IQR}} and
#' \eqn{Q(0.75) + k \cdot \mathrm{IQR}} (Tukey, 1977).
#'
#' Missing values are not removed: an `NA` or `NaN` anywhere makes every
#' result `NaN` (where [stats::quantile()] stops). An empty `x` gives
#' `NaN`. `cqnt()` refuses an empty `probs` and any probability outside
#' \eqn{[0, 1]} or missing.
#'
#' @param x Numeric (or logical) vector.
#' @param probs Probabilities, each in \eqn{[0, 1]}.
#' @param const The scale factor of `cmad()`, one number from 0: the
#'   default 1.4826 is the rounded constant [stats::mad()] uses, so the two
#'   agree exactly; `1` gives the unscaled median deviation.
#' @param k The fence width in IQRs, one number from 0 (1.5 is Tukey's
#'   convention, 3 the usual far-out cut).
#' @return `cqnt()`: a numeric vector the length of `probs`, named as
#'   [stats::quantile()] names it. `cmed()`, `cmad()`, `ciqr()`: a number.
#'   `ctkyfnc()`: a named numeric vector, `lower` and `upper`.
#' @references Hyndman RJ, Fan Y (1996). Sample quantiles in statistical
#'   packages. *The American Statistician* 50(4), 361--365.
#'   \doi{10.1080/00031305.1996.10473566}
#'
#' Tukey JW (1977). *Exploratory Data Analysis*. Addison-Wesley.
#' @seealso [cmom()], [ctrmmean()].
#' @examples
#' x <- c(2, 4, 4, 4, 5, 5, 7, 9)
#' cqnt(x, c(0.25, 0.5, 0.75))
#' cmed(x)
#' cmad(x)
#' all.equal(cmad(x), stats::mad(x))
#' ciqr(x)
#' ctkyfnc(x)
#' wild <- c(x, 1000)
#' f <- ctkyfnc(wild)
#' wild[wild < f[["lower"]] | wild > f[["upper"]]]
#' cmed(c(1, NA, 3))
#' try(cqnt(x, 1.5))
#' @name crob
#' @export
cqnt <- function(x, probs = c(0, 0.25, 0.5, 0.75, 1)) {
  x <- .pnum(x, "x", "cqnt")
  probs <- .pnum(probs, "probs", "cqnt")
  if (!length(probs)) {
    stop("cqnt: `probs` is empty.", call. = FALSE)
  }
  if (anyNA(probs) || any(probs < 0 | probs > 1)) {
    stop("cqnt: `probs` must lie in [0, 1].", call. = FALSE)
  }
  out <- .pcall(6L, x, probs)
  names(out) <- paste0(formatC(100 * probs, format = "fg", width = 1,
                               digits = 7L), "%")
  out
}

#' @rdname crob
#' @export
cmed <- function(x) .pcall(6L, .pnum(x, "x", "cmed"), 0.5)

#' @rdname crob
#' @export
cmad <- function(x, const = 1.4826) {
  x <- .pnum(x, "x", "cmad")
  .pcall(7L, x, a = .pone(const, "const", "cmad"))
}

#' @rdname crob
#' @export
ciqr <- function(x) {
  q <- .pcall(6L, .pnum(x, "x", "ciqr"), c(0.25, 0.75))
  q[2L] - q[1L]
}

#' @rdname crob
#' @export
ctkyfnc <- function(x, k = 1.5) {
  x <- .pnum(x, "x", "ctkyfnc")
  k <- .pone(k, "k", "ctkyfnc")
  q <- .pcall(6L, x, c(0.25, 0.75))
  c(lower = q[1L] - k * (q[2L] - q[1L]), upper = q[2L] + k * (q[2L] - q[1L]))
}

#' Trimmed and winsorized means (C backend)
#'
#' Two ways to stop a handful of extreme rows dominating a column's
#' centre. `ctrmmean()` DISCARDS the `floor(n * trim)` largest and
#' smallest values, as `mean(x, trim = )` does. `cwinsmean()` instead
#' PULLS THEM IN to the most extreme surviving values, so every
#' observation still counts, usually the better choice when the extremes
#' are real measurements rather than errors.
#'
#' With \eqn{g = \lfloor n \cdot \mathrm{trim} \rfloor} and the sorted
#' values \eqn{x_{(1)} \le \dots \le x_{(n)}}, the trimmed mean is
#' \eqn{\frac{1}{n - 2g} \sum_{i=g+1}^{n-g} x_{(i)}} and the winsorized
#' mean replaces \eqn{x_{(1)}, \dots, x_{(g)}} by \eqn{x_{(g+1)}} and
#' \eqn{x_{(n-g+1)}, \dots, x_{(n)}} by \eqn{x_{(n-g)}} before averaging all
#' \eqn{n} values. At `trim = 0.5` both return the median, as `mean()`
#' does. An empty `x`, or any `NA` or `NaN`, gives `NaN`.
#'
#' @param x Numeric (or logical) vector.
#' @param trim The proportion trimmed (or winsorized) at each end, one
#'   number from 0 to 0.5.
#' @return A number in the units of `x`.
#' @seealso [cmed()], [cmad()], [cmom()].
#' @examples
#' x <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 100)
#' mean(x)
#' ctrmmean(x, 0.1)
#' cwinsmean(x, 0.1)
#' all.equal(ctrmmean(x, 0.2), mean(x, trim = 0.2))
#' all.equal(ctrmmean(x, 0.5), cmed(x))
#' try(cwinsmean(x, 0.7))
#' @name ctrm
#' @export
ctrmmean <- function(x, trim = 0.1) {
  x <- .pnum(x, "x", "ctrmmean")
  .pcall(8L, x, a = .pone(trim, "trim", "ctrmmean", 0, 0.5))
}

#' @rdname ctrm
#' @export
cwinsmean <- function(x, trim = 0.1) {
  x <- .pnum(x, "x", "cwinsmean")
  .pcall(9L, x, a = .pone(trim, "trim", "cwinsmean", 0, 0.5))
}

#' Weighted mean and variance (C backend)
#'
#' The weighted mean and variance of a numeric vector, for data where some
#' observations count for more than others: an average over areas
#' weighted by population, or measurements weighted by their precision.
#'
#' The mean is \eqn{\bar x_w = \sum_i w_i x_i / \sum_i w_i}. The weights
#' are read as RELIABILITY weights, so the variance is
#' \deqn{s_w^2 = \frac{\sum_i w_i (x_i - \bar x_w)^2}
#'   {\sum_i w_i - \sum_i w_i^2 / \sum_i w_i},}
#' the unbiased estimator under that reading; with every weight 1 it is
#' [stats::var()]. Frequency weights (a weight of 3 meaning three
#' identical rows) call for \eqn{\sum w - 1} instead, and survey design
#' weights for a design-based variance: do not use this variance for
#' those. An `NA` in either vector gives `NaN` for both entries.
#'
#' @param x Numeric vector of observations.
#' @param w Non-negative finite weights, the length of `x`, not all zero.
#' @return A named numeric vector: `mean` and `variance`.
#' @seealso [cmom()], [cipw()].
#' @examples
#' x <- c(10, 20, 30, 40)
#' w <- c(1, 1, 2, 4)
#' cwtd(x, w)
#' all.equal(cwtd(x, w)[["mean"]], stats::weighted.mean(x, w))
#' all.equal(cwtd(x, rep(1, 4))[["variance"]], stats::var(x))
#' try(cwtd(x, c(1, -1, 1, 1)))
#' @export
cwtd <- function(x, w) {
  x <- .pnum(x, "x", "cwtd")
  w <- .pnum(w, "w", "cwtd")
  .psame(x, w, "x", "w", "cwtd")
  if (any(w < 0, na.rm = TRUE) || any(is.infinite(w)) ||
        !any(w > 0, na.rm = TRUE)) {
    stop("cwtd: `w` must be non-negative and finite, not all zero.",
         call. = FALSE)
  }
  out <- .pcall(10L, x, w)
  names(out) <- c("mean", "variance")
  out
}

#' Rank correlation and midranks (C backend)
#'
#' `cspr()` is Spearman's rho, the Pearson correlation of the ranks: it
#' measures monotone association and is unaffected by any order-preserving
#' transformation of either variable. `cmidrk()` gives the ranks
#' themselves; tied values share the average of the ranks they span, which
#' is what makes the result agree with [stats::cor()] on tied data.
#'
#' With \eqn{r_i} and \eqn{s_i} the midranks of \eqn{x_i} and \eqn{y_i},
#' \deqn{\rho = \frac{\sum_i (r_i - \bar r)(s_i - \bar s)}
#'   {\sqrt{\sum_i (r_i - \bar r)^2 \sum_i (s_i - \bar s)^2}}}
#' (Spearman, 1904); the shortcut \eqn{1 - 6 \sum d_i^2 / (n(n^2 - 1))}
#' is not used, so ties are handled correctly. `cspr()` gives `NaN` for an
#' `NA` in either vector and for a constant vector. `cmidrk()` keeps an
#' `NA` in place and ranks the rest, as `rank(x, na.last = "keep")` does.
#'
#' @param x,y Numeric (or logical) vectors of the same length.
#' @return `cspr()`: a number in \eqn{[-1, 1]}. `cmidrk()`: the rank of
#'   each value, the length of `x`.
#' @references Spearman C (1904). The proof and measurement of association
#'   between two things. *The American Journal of Psychology* 15(1),
#'   72--101. \doi{10.2307/1412159}
#' @seealso [ccorr()].
#' @examples
#' x <- c(1, 2, 3, 4, 5)
#' y <- c(2, 4, 9, 16, 25)
#' cspr(x, y)
#' ccorr(x, y)
#' xt <- c(1, 2, 2, 2, 5, 5, 7)
#' yt <- c(3, 1, 1, 4, 4, 9, 2)
#' all.equal(cspr(xt, yt), stats::cor(xt, yt, method = "spearman"))
#' cmidrk(xt)
#' cmidrk(c(3, NA, 1, 3))
#' try(cspr(1:3, 1:4))
#' @name crank
#' @export
cspr <- function(x, y) {
  x <- .pnum(x, "x", "cspr")
  y <- .pnum(y, "y", "cspr")
  .psame(x, y, "x", "y", "cspr")
  .pcall(11L, x, y)
}

#' @rdname crank
#' @export
cmidrk <- function(x) {
  x <- .pnum(x, "x", "cmidrk")
  out <- rep(NA_real_, length(x))
  ok <- !is.na(x)
  out[ok] <- .pcall(12L, x[ok])
  out
}

#' Covariance matrix (C backend)
#'
#' The covariance of every pair of columns of a numeric matrix or data
#' frame. Column names are carried to both dimensions of the result.
#'
#' Entry \eqn{(j, k)} is \eqn{\frac{1}{n - 1} \sum_i (x_{ij} - \bar
#' x_j)(x_{ik} - \bar x_k)}, as [stats::cov()] gives it; the diagonal holds
#' the column variances. Missing values are not removed: an `NA` in a
#' column makes every covariance with that column `NA`. A single row gives
#' `NaN` throughout. A plain vector, and a data frame with a column that
#' is not numeric, are refused.
#'
#' @param x A numeric matrix, or a data frame of numeric columns (rows
#'   observations, columns variables).
#' @return A symmetric `ncol(x)` by `ncol(x)` matrix.
#' @seealso [ccorr()], [cmom()].
#' @examples
#' X <- cbind(a = c(1, 2, 3, 4), b = c(2, 4, 7, 8), c = c(5, 3, 2, 1))
#' ccovm(X)
#' all.equal(ccovm(X), stats::cov(X))
#' ccovm(data.frame(u = 1:5, v = c(2, 1, 4, 3, 6)))
#' try(ccovm(1:5))
#' @export
ccovm <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, NA))) {
      stop("ccovm: every column of `x` must be numeric.", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("ccovm: `x` must be a numeric matrix or data frame.", call. = FALSE)
  }
  out <- matrix(.pcall(13L, as.numeric(x), a = ncol(x)), ncol(x))
  if (!is.null(colnames(x))) {
    dimnames(out) <- list(colnames(x), colnames(x))
  }
  out
}

#' Bootstrap replicate means (C backend)
#'
#' `B` resamples of `x`, drawn with replacement and each the length of
#' `x`, with the mean of each. The spread of these replicate means
#' estimates the standard error of the sample mean, and their quantiles
#' give a percentile interval, without assuming the data are normal (Efron
#' and Tibshirani, 1993). The draws come from the core's own 64-bit
#' Mersenne Twister seeded by `seed`, not from R's generator, so a seed
#' gives the same replicates in every front end and on every platform, and
#' R's random stream is left untouched.
#'
#' Replicate \eqn{b} draws \eqn{n} indices uniformly from \eqn{1, \dots,
#' n} (`std::mt19937_64` seeded with `seed`, an unbiased draw by Lemire's
#' multiply-and-reject) and returns the mean of the values drawn,
#' accumulated in extended precision. The replicates are one stream: the
#' first `B` for a seed are the same whatever larger `B` is asked for. An
#' empty `x` gives `B` values of `NaN`; an `NA` in `x` makes every
#' replicate that draws it `NA`. A long run can be interrupted.
#'
#' @param x Numeric (or logical) vector to resample.
#' @param B The number of replicates, one whole number from 1 to
#'   67108864 (2^26).
#' @param seed The seed of the core's generator, one whole number from 0,
#'   below \eqn{2^{53}}.
#' @return A numeric vector of `B` replicate means.
#' @references Efron B, Tibshirani RJ (1993). *An Introduction to the
#'   Bootstrap*. Chapman and Hall.
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(50, mean = 5)
#' reps <- cbootmean(x, B = 500, seed = 7)
#' c(sample = mean(x), bootstrap = mean(reps))
#' stats::quantile(reps, c(0.025, 0.975))
#' identical(cbootmean(x, 10, seed = 1), cbootmean(x, 100, seed = 1)[1:10])
#' try(cbootmean(x, B = 0))
#' @export
cbootmean <- function(x, B = 1000L, seed = 42L) {
  x <- .pnum(x, "x", "cbootmean")
  B <- .pone(B, "B", "cbootmean", 1, .Machine$integer.max, TRUE)
  .pmax(B, "B", "cbootmean")
  .pcall(14L, x, a = c(B, .pseed(seed, "cbootmean")))
}

#' Trimmed inverse-probability weights (C backend)
#'
#' The inverse-probability-of-treatment weights \eqn{1/e} for the treated
#' and \eqn{1/(1-e)} for the untreated, with the propensity score clamped
#' into `[trimlo, trimhi]` first: an untrimmed score near 0 or 1 gives one
#' observation a weight large enough to dominate the estimate. The
#' weights reweight a treated and an untreated group so that each
#' resembles the whole sample on the covariates behind the propensity
#' score (Rosenbaum and Rubin, 1983).
#'
#' With the clamped score \eqn{\tilde e_i = \min(\max(e_i, l), u)}, the
#' weight is \eqn{1 / \tilde e_i} when \eqn{T_i = 1} and \eqn{1 / (1 -
#' \tilde e_i)} when \eqn{T_i = 0}; with the defaults none exceeds 100. A
#' treatment coded other than 0/1, a propensity outside \eqn{[0, 1]} and
#' bounds without \eqn{0 < l < u < 1} are refused. A missing treatment or
#' propensity gives a missing weight for that unit only. The weights are
#' not normalised.
#'
#' @param treat Treatment indicator, 0/1 or `TRUE`/`FALSE`; `NA` allowed.
#' @param ps Propensity scores in \eqn{[0, 1]}, the length of `treat`;
#'   `NA` allowed.
#' @param trimlo,trimhi The clamp bounds, `0 < trimlo < trimhi < 1`.
#' @return The weights, the length of `treat`.
#' @references Rosenbaum PR, Rubin DB (1983). The central role of the
#'   propensity score in observational studies for causal effects.
#'   *Biometrika* 70(1), 41--55. \doi{10.1093/biomet/70.1.41}
#' @seealso [cwtd()].
#' @examples
#' cipw(c(1, 0, 1, 0), c(0.5, 0.25, 0.02, 0.9))
#' cipw(c(1, 0, 1, 0), c(0.5, 0.25, 0.02, 0.9), trimlo = 0.1)
#' cipw(c(TRUE, FALSE), c(0.4, 0.4))
#' try(cipw(c(1, 2), c(0.5, 0.5)))
#' @export
cipw <- function(treat, ps, trimlo = 0.01, trimhi = 0.99) {
  treat <- .pnum(treat, "treat", "cipw")
  ps <- .pnum(ps, "ps", "cipw")
  .psame(treat, ps, "treat", "ps", "cipw")
  lo <- .pone(trimlo, "trimlo", "cipw", 0, 1)
  hi <- .pone(trimhi, "trimhi", "cipw", 0, 1)
  if (!(lo > 0 && hi < 1 && lo < hi)) {
    stop("cipw: need 0 < `trimlo` < `trimhi` < 1.", call. = FALSE)
  }
  if (!all(treat[!is.na(treat)] %in% c(0, 1))) {
    stop("cipw: `treat` must be coded 0/1.", call. = FALSE)
  }
  if (any(ps < 0 | ps > 1, na.rm = TRUE)) {
    stop("cipw: `ps` must lie in [0, 1].", call. = FALSE)
  }
  out <- rep(NA_real_, length(treat))
  k <- !is.na(treat) & !is.na(ps)
  out[k] <- .pcall(15L, treat[k], ps[k], c(lo, hi))
  out
}

#' Regularized incomplete gamma function (C backend)
#'
#' The lower regularized incomplete gamma function \eqn{P(a, x)}, the
#' distribution function of a gamma distribution with shape `a` and unit
#' rate; the p-value of a chi-square statistic \eqn{q} on \eqn{k} degrees
#' of freedom is \eqn{1 - P(k/2, q/2)}.
#'
#' \deqn{P(a, x) = \frac{1}{\Gamma(a)} \int_0^x t^{a-1} e^{-t} \, dt,}
#' computed by the package's own gamma kernel, which agrees with
#' `stats::pgamma(x, shape)` to rounding. For a rate \eqn{r}, pass `x * r`.
#' A negative `x` gives 0, `NA` gives `NA` and `Inf` gives 1; `shape = Inf`
#' gives 0 for every finite `x`, as `pgamma()` does.
#'
#' @param shape The shape \eqn{a}, one positive number.
#' @param x Numeric vector of quantiles.
#' @return A numeric vector of probabilities, the length of `x`.
#' @seealso [stats::pgamma()].
#' @examples
#' cgamcdf(3.5, c(0.5, 1, 4, 12))
#' all.equal(cgamcdf(3.5, c(0.5, 1, 4, 12)),
#'           stats::pgamma(c(0.5, 1, 4, 12), shape = 3.5))
#' 1 - cgamcdf(2 / 2, 5.99 / 2)
#' try(cgamcdf(0, 1))
#' @export
cgamcdf <- function(shape, x) {
  if (!is.numeric(shape) || length(shape) != 1L || is.na(shape) ||
        shape <= 0) {
    stop("cgamcdf: `shape` must be one positive number.", call. = FALSE)
  }
  x <- .pnum(x, "x", "cgamcdf")
  out <- .pcall(17L, x, a = shape)
  out[is.na(x)] <- x[is.na(x)]
  out
}

#' Reproducible uniforms (splitmix64)
#'
#' Returns `n` pseudo-random numbers on [0, 1) from a 64-bit seed with
#' the splitmix64 generator. The stream depends only on `seed`, and any
#' program that implements splitmix64 (Python, Go, a spreadsheet macro)
#' computes the same numbers, so a subsample or a jitter drawn in one
#' language can be reproduced exactly in another.
#'
#' splitmix64 (Steele, Lea and Flood, 2014) advances a 64-bit state by the
#' constant `0x9E3779B97F4A7C15` and mixes it with two multiply-xorshift
#' rounds; the top 53 bits of each output divided by \eqn{2^{53}} give a
#' double in [0, 1). R's own generator and `set.seed()` are neither used
#' nor changed. It is not a statistical replacement for R's generators.
#'
#' @param n The number of uniforms, one whole number from 0 to
#'   67108864 (2^26).
#' @param seed The seed, one whole number from 0, below \eqn{2^{53}} (so
#'   that it is exact as a double).
#' @return A numeric vector of `n` values in [0, 1).
#' @references Steele GL, Lea D, Flood CH (2014). Fast splittable
#'   pseudorandom number generators. *OOPSLA 2014*, 453--472.
#'   \doi{10.1145/2660193.2660195}
#' @examples
#' cunif(3, 42)
#' all.equal(cunif(3, 7), cunif(10, 7)[1:3])
#' set.seed(1); a <- runif(1)
#' set.seed(1); invisible(cunif(100, 9)); b <- runif(1)
#' identical(a, b)
#' try(cunif(-1, 1))
#' try(cunif(2, 2^53))
#' @export
cunif <- function(n, seed) {
  n <- .pone(n, "n", "cunif", 0, .Machine$integer.max, TRUE)
  .pmax(n, "n", "cunif")
  .pcall(16L, numeric(), a = c(n, .pseed(seed, "cunif")))
}

.pcall <- function(op, x, y = NULL, a = NULL) {
  .Call(C_clab_r_prim, op, x, y, if (!is.null(a)) as.numeric(a))
}

# numbers or TRUE and FALSE, as doubles; NULL is no values
.pnum <- function(x, nm, fn) {
  if (is.null(x)) {
    return(numeric())
  }
  if (!is.atomic(x) || !(is.numeric(x) || is.logical(x)) || is.factor(x)) {
    stop(fn, ": `", nm, "` must be numbers (or TRUE and FALSE).",
         call. = FALSE)
  }
  as.numeric(x)
}

# one number from lo (to hi), finite
.pone <- function(x, nm, fn, lo = 0, hi = Inf, whole = FALSE) {
  ok <- is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lo &&
    x <= hi && (!whole || x == floor(x))
  if (!ok) {
    what <- if (whole) "one whole number" else "one number"
    to <- if (is.finite(hi) && !whole) paste0(" to ", format(hi)) else ""
    stop(fn, ": `", nm, "` must be ", what, " from ", format(lo), to, ".",
         call. = FALSE)
  }
  as.numeric(x)
}

# A count that sizes the result, at most 2^26 (the bound of the kernels).
.pmax <- function(x, nm, fn) {
  if (x > 2^26) {
    stop(fn, ": `", nm, "` must be at most 67108864.", call. = FALSE)
  }
  invisible(x)
}

.psame <- function(x, y, a, b, fn) {
  if (length(x) != length(y)) {
    stop(fn, ": `", a, "` and `", b, "` must have the same length.",
         call. = FALSE)
  }
}

.pseed <- function(seed, fn) {
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || seed < 0 ||
        seed != floor(seed) || seed >= 2^53) {
    stop(fn, ": `seed` must be one whole number from 0, below 2^53.",
         call. = FALSE)
  }
  as.numeric(seed)
}

.pnorm <- function(x, mean, sd, lg, fn) {
  x <- .pnum(x, "x", fn)
  if (!is.numeric(mean) || length(mean) != 1L || !is.finite(mean)) {
    stop(fn, ": `mean` must be one finite number.", call. = FALSE)
  }
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0) {
    stop(fn, ": `sd` must be one positive finite number.", call. = FALSE)
  }
  .pcall(4L, x, a = c(mean, sd, lg))
}
