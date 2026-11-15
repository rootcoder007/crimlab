# Scaling and ideal points: thin wrappers over the kernels of src/svt.cpp
# (src/svk.h). Each passes its arguments to the kernel table in
# src/svla.cpp, whose defaults the wrappers repeat, and classes the result.

.sv <- function(fn, args, cls = "svfit") {
  args <- lapply(args, function(a) if (is.data.frame(a)) as.matrix(a) else a)
  r <- .Call(C_clab_r_svk, fn, args)
  class(r) <- c(fn, cls)
  r
}

#' Vote agreement between legislators
#'
#' The share of the roll calls on which each pair of legislators (judges,
#' council members) voted the same way, over the roll calls both took
#' part in (Armstrong et al. 2021, section 3.2.2): the similarity matrix
#' that metric scaling of a chamber starts from.
#'
#' @param x A legislator-by-vote matrix of votes (any coding; equal
#'   values agree), `NA` for an absence. A vector is one vote.
#' @return A list of class `votagr`: `agree` (legislators by legislators,
#'   1 on the diagonal, `NA` for a pair with no vote in common), `mean`
#'   (the mean agreement over the pairs), `n` and `m` (legislators and
#'   votes).
#' @references Armstrong, D. A., Bakker, R., Carroll, R., Hare, C.,
#'   Poole, K. T. and Rosenthal, H. (2021). *Analyzing Spatial Models of
#'   Choice and Judgment*, 2nd ed. Chapman & Hall/CRC.
#' @examples
#' # five judges on eight simulated cases
#' set.seed(1)
#' v <- matrix(rbinom(40, 1, 0.5), 5, 8)
#' v[2, 3] <- NA
#' votagr(v)$agree
#' @family scaling
#' @export
votagr <- function(x) .sv("votagr", list(x = x))

#' Aldrich-McKelvey scaling
#'
#' Recovers the positions of stimuli (parties, candidates, courts) on one
#' scale from respondents' placements of them, each respondent allowed
#' an own shift and stretch of the scale: \eqn{z_{ij} = a_i + b_i \hat
#' z_j + \epsilon_{ij}}. The positions are the least-squares solution of
#' Aldrich and McKelvey (1977): with \eqn{X_i = [1, z_i]} and \eqn{P_i}
#' the projection onto its columns, \eqn{\hat z} is the eigenvector of
#' \eqn{\sum_i (I - P_i)} with the smallest eigenvalue orthogonal to the
#' constant. Respondents with a missing placement, or who place every
#' stimulus at the same point, are left out of that sum. The positions
#' are standardised to mean 0 and standard deviation 1, the first
#' stimulus on the left. The respondent intercepts and slopes are each
#' respondent's least-squares regression of the placements on those
#' positions (defined with two placements or more; a slope of zero is
#' set to 1e-10).
#'
#' @param Z A respondent-by-stimulus matrix of placements, `NA` for a
#'   missing one.
#' @return A list of class `amscal`: `zhat` (stimulus positions),
#'   `alpha` and `beta` (each respondent's intercept and slope) and
#'   `weights` (\eqn{|b_i|} scaled to sum to the number of respondents).
#' @references Aldrich, J. H. and McKelvey, R. D. (1977). A method of
#'   scaling with applications to the 1968 and 1972 presidential
#'   elections. *American Political Science Review*, 71(1), 111-130.
#'   \doi{10.2307/1956957}
#'
#'   Poole, K. T. (1998). Recovering a basic space from a set of issue
#'   scales. *American Journal of Political Science*, 42(3), 954-993.
#'   \doi{10.2307/2991737}
#' @examples
#' # 40 simulated respondents place five stimuli on a 1-7 scale
#' set.seed(1)
#' truth <- c(-1.2, -0.4, 0.1, 0.9, 1.5)
#' Z <- t(sapply(1:40, function(i) 4 + rnorm(1, 0, 0.3) +
#'   runif(1, 0.6, 1.4) * truth + rnorm(5, 0, 0.3)))
#' round(amscal(Z)$zhat, 2)
#' round((truth - mean(truth)) / sd(truth), 2)
#' @family scaling
#' @export
amscal <- function(Z) .sv("amscal", list(Z = Z))

#' Basic Space (blackbox) scaling
#'
#' Recovers respondents' positions from a respondent-by-issue matrix of
#' responses by the decomposition of Poole (1998), \eqn{X_0 = \Psi W' +
#' J_n c' + E_0} over the observed cells. Each dimension starts from the
#' sign-corrected respondent means (column polarities from the pairwise
#' correlations) and is refined by four alternating least-squares sweeps
#' on the residuals of the earlier dimensions; up to five joint sweeps
#' follow, stopping once the two half-step error sums of squares differ
#' by less than 0.01; the fitted \eqn{\Psi W'} is then re-expressed by
#' its singular value decomposition, \eqn{\Psi = U D^{1/2}}, \eqn{W = V
#' D^{1/2}}. Respondents with fewer than `minscale` responses are not
#' scaled. The sign of each dimension is arbitrary.
#'
#' @param X A respondent-by-issue matrix of responses, `NA` for missing.
#' @param ndim The number of dimensions (at most one fewer than the
#'   issues and than the scaled respondents).
#' @param minscale The fewest responses a respondent needs to be scaled
#'   (capped at the number of issues).
#' @return A list of class `bbscal`: `ideal` (respondents by dimensions,
#'   `NA` rows for those not scaled), `weights` (issues by dimensions),
#'   `eigen` and `singular` (of the fitted \eqn{\Psi W'}), `fit` (\eqn{1
#'   - SSE/SS}, SS about the grand mean of the observed cells), `cmeans`
#'   (the issue intercepts \eqn{c}) and `ndim`.
#' @references Poole, K. T. (1998). Recovering a basic space from a set of
#'   issue scales. *American Journal of Political Science*, 42(3),
#'   954-993. \doi{10.2307/2991737}
#' @examples
#' # 30 simulated respondents on eight issues, ten answers missing
#' set.seed(1)
#' X <- matrix(rnorm(60), 30, 2) %*% matrix(rnorm(16), 2, 8) +
#'   matrix(rnorm(240, 0, 0.3), 30, 8)
#' X[sample(240, 10)] <- NA
#' bbscal(X, ndim = 2)$fit
#' @family scaling
#' @export
bbscal <- function(X, ndim = 2, minscale = 8) {
  .sv("bbscal", list(X = X, ndim = ndim, minscale = minscale))
}

#' Optimal Classification scaling
#'
#' Nonparametric ideal points from binary votes in the spirit of Poole
#' (2000): each vote's cutting direction is the difference of the mean
#' positions of its yea and nay voters, each legislator moves to the
#' unit-length sum of the directions of the sides taken, for `maxit`
#' sweeps from each of `nstart` random starts; the start with the fewest
#' classification errors (a vote is split at the midpoint of the two
#' sides' mean projections) is kept.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @param ndim The number of dimensions.
#' @param maxit The sweeps per start.
#' @param nstart The random starts.
#' @param seed The seed of the starts, as [set.seed()] takes it.
#' @return A list of class `optcls`: `ideal`, `normals` (votes by
#'   dimensions), `pre` (the proportional reduction in error), `errors`,
#'   `null` (the errors of voting with the majority on every vote) and
#'   `ndim`.
#' @references Poole, K. T. (2000). Nonparametric unfolding of binary
#'   choice data. *Political Analysis*, 8(3), 211-237.
#'   \doi{10.1093/oxfordjournals.pan.a029814}
#' @examples
#' # 20 simulated legislators, 30 votes cut along one dimension
#' set.seed(1)
#' x <- sort(rnorm(20))
#' v <- outer(x, rnorm(30), ">") * 1
#' optcls(v, maxit = 20, nstart = 2)$pre
#' @family scaling
#' @export
optcls <- function(votes, ndim = 1, maxit = 500, nstart = 10, seed = 42) {
  .sv(
    "optcls",
    list(votes = votes, ndim = ndim, maxit = maxit, nstart = nstart,
         seed = seed)
  )
}

#' Double-centre a distance matrix
#'
#' \eqn{B = -\frac{1}{2} H D^{(2)} H} with \eqn{H = I - n^{-1} 1 1'} and
#' \eqn{D^{(2)}} the squared distances: the inner-product matrix of
#' classical scaling (Torgerson 1952).
#'
#' @param D A symmetric matrix of distances.
#' @return A list of class `dblctr` with the matrix `B`.
#' @references Torgerson, W. S. (1952). Multidimensional scaling: I.
#'   Theory and method. *Psychometrika*, 17(4), 401-419.
#'   \doi{10.1007/BF02288916}
#' @examples
#' D <- as.matrix(dist(cbind(c(0, 1, 3, 4), c(0, 1, 0, 2))))
#' dblctr(D)$B
#' @family scaling
#' @export
dblctr <- function(D) .sv("dblctr", list(D = D))

#' Classical (metric) multidimensional scaling
#'
#' Torgerson scaling: the leading eigenvectors of the double-centred
#' squared distances ([dblctr()]), each scaled by the square root of its
#' eigenvalue (negative ones set to 0). The fit is Mardia's criterion,
#' the kept positive eigenvalues over the sum of all absolute ones; the
#' stress is Kruskal's stress-1 against the positive distances. The sign
#' of each dimension is arbitrary.
#'
#' @param D A symmetric matrix of distances.
#' @param ndim The number of dimensions.
#' @return A list of class `cmdsf`: `coords`, `eigen` (the first `ndim`
#'   eigenvalues), `stress`, `fit` and `B`.
#' @references Torgerson, W. S. (1952). Multidimensional scaling: I.
#'   Theory and method. *Psychometrika*, 17(4), 401-419.
#'   \doi{10.1007/BF02288916}
#'
#'   Mardia, K. V. (1978). Some properties of classical
#'   multi-dimensional scaling. *Communications in Statistics - Theory and
#'   Methods*, 7(13), 1233-1241.
#' @examples
#' set.seed(1)
#' D <- as.matrix(dist(matrix(rnorm(40), 10)))
#' cmdsf(D)$fit
#' @family scaling
#' @export
cmdsf <- function(D, ndim = 2) .sv("cmdsf", list(D = D, ndim = ndim))

#' SMACOF stress minimisation
#'
#' Metric scaling by majorisation (de Leeuw 1977): each step is the
#' Guttman transform \eqn{X \leftarrow V^{+} B(X) X} with \eqn{V =
#' \mathrm{diag}(W 1) - W}, until the raw stress \eqn{\frac{1}{2} \sum
#' w_{ij} (d_{ij} - \delta_{ij})^2} changes by less than `tol`.
#'
#' @param D A symmetric matrix of dissimilarities.
#' @param ndim The number of dimensions.
#' @param maxit The most steps.
#' @param tol The convergence tolerance on the change in stress.
#' @param w `NULL` for equal weights, or a matrix of weights.
#' @param init `NULL` for a random start (normal, from seed 42 as
#'   [set.seed()] takes it), or a starting configuration.
#' @return A list of class `smacf`: `coords`, `stress`, `iter` and
#'   `conv`.
#' @references de Leeuw, J. (1977). Applications of convex analysis to
#'   multidimensional scaling. In J. R. Barra et al. (eds.), *Recent
#'   Developments in Statistics*, 133-145. North-Holland.
#'
#'   Borg, I. and Groenen, P. J. F. (2005). *Modern Multidimensional
#'   Scaling*, 2nd ed. Springer. \doi{10.1007/0-387-28981-X}
#' @examples
#' set.seed(1)
#' D <- as.matrix(dist(matrix(rnorm(40), 10)))
#' smacf(D)$stress
#' @family scaling
#' @export
smacf <- function(D, ndim = 2, maxit = 300, tol = 1e-6, w = NULL,
                  init = NULL) {
  .sv(
    "smacf",
    list(D = D, ndim = ndim, maxit = maxit, tol = tol, w = w, init = init)
  )
}

#' Nonmetric multidimensional scaling
#'
#' Kruskal's nonmetric scaling: at each step the distances of the current
#' configuration are fitted by a monotone (isotonic, pool-adjacent-
#' violators) regression on the order of the dissimilarities, and one
#' SMACOF step ([smacf()]) is taken towards the fitted disparities, until
#' no coordinate moves by `tol`. The stress is Kruskal's stress-1 of the
#' last step. The start is normal from seed 42, as [set.seed()] takes it.
#'
#' @inheritParams smacf
#' @param tol The convergence tolerance on the largest change of a
#'   coordinate.
#' @return A list of class `nmdsf`: `coords`, `stress`, `iter` and
#'   `conv`.
#' @references Kruskal, J. B. (1964). Nonmetric multidimensional scaling:
#'   a numerical method. *Psychometrika*, 29(2), 115-129.
#'   \doi{10.1007/BF02289694}
#' @examples
#' set.seed(1)
#' D <- as.matrix(dist(matrix(rnorm(40), 10)))
#' nmdsf(D, maxit = 50)$stress
#' @family scaling
#' @export
nmdsf <- function(D, ndim = 2, maxit = 300, tol = 1e-6) {
  .sv("nmdsf", list(D = D, ndim = ndim, maxit = maxit, tol = tol))
}

#' Fit of a scaling by its eigenvalues
#'
#' Mardia's criterion: the share of each dimension, its eigenvalue (0 when
#' negative) over the sum of the absolute values of all eigenvalues, and
#' the cumulative shares.
#'
#' @param ev The eigenvalues.
#' @return A list of class `mdsfit`: `fit`, `cumfit` and `ev` (empty fits
#'   when every eigenvalue is 0).
#' @references Mardia, K. V. (1978). Some properties of classical
#'   multi-dimensional scaling. *Communications in Statistics - Theory and
#'   Methods*, 7(13), 1233-1241.
#' @examples
#' mdsfit(c(4, 2, 1, 0.5))
#' @family scaling
#' @export
mdsfit <- function(ev) .sv("mdsfit", list(ev = ev))

#' Stress of an unfolding
#'
#' The weighted sum of squared differences between the distances from
#' respondents to stimuli in a configuration and the observed
#' dissimilarities, over the observed cells.
#'
#' @param xr Respondent coordinates (respondents by dimensions).
#' @param xs Stimulus coordinates (stimuli by dimensions).
#' @param D The respondent-by-stimulus dissimilarities, `NA` for missing.
#' @param w `NULL` for equal weights, or a matrix of weights like `D`.
#' @return A list of class `unfstr` with `stress`.
#' @references Coombs, C. H. (1964). *A Theory of Data*. Wiley.
#' @examples
#' set.seed(1)
#' unfstr(matrix(rnorm(6), 3), matrix(rnorm(8), 4), matrix(runif(12), 3))
#' @family scaling
#' @export
unfstr <- function(xr, xs, D, w = NULL) {
  .sv("unfstr", list(xr = xr, xs = xs, D = D, w = w))
}

#' Least-squares metric unfolding
#'
#' Unfolding of a respondent-by-stimulus matrix in the alternating
#' least-squares form of MLSMU6 (Poole 1984): respondents and stimuli
#' take gradient steps in turn against the row-centred dissimilarities,
#' each step scaled by the largest eigenvalue of the gradient's
#' cross-product and both clouds centred, until the stress
#' ([unfstr()]) changes by a relative `tol`; of `nstart` normal starts
#' (seed 42, as [set.seed()] takes it) the lowest stress is kept.
#'
#' @param D A respondent-by-stimulus matrix of dissimilarities.
#' @param ndim The number of dimensions.
#' @param maxit The most sweeps per start.
#' @param tol The convergence tolerance on the relative change in stress.
#' @param nstart The random starts.
#' @return A list of class `mlsmu`: `xr`, `xs`, `stress`, `iter` and
#'   `conv`.
#' @references Poole, K. T. (1984). Least squares metric, unidimensional
#'   unfolding. *Psychometrika*, 49(3), 311-323.
#'
#'   Bakker, R. and Poole, K. T. (2013). Bayesian metric
#'   multidimensional scaling. *Political Analysis*, 21(1), 125-140.
#'   \doi{10.1093/pan/mps039}
#' @examples
#' set.seed(1)
#' D <- matrix(runif(120), 20, 6)
#' mlsmu(D, nstart = 1, maxit = 50)$stress
#' @family scaling
#' @export
mlsmu <- function(D, ndim = 2, maxit = 200, tol = 1e-6, nstart = 5) {
  .sv(
    "mlsmu",
    list(D = D, ndim = ndim, maxit = maxit, tol = tol, nstart = nstart)
  )
}

#' SMACOF unfolding
#'
#' Unfolding by majorisation: respondents and stimuli are embedded in one
#' space by SMACOF steps ([smacf()]) on the joint configuration with
#' weight 1 between a respondent and a stimulus and 0 within each set
#' (Borg and Groenen 2005, chapter 14). The start is normal from seed 42,
#' as [set.seed()] takes it.
#'
#' @param D A respondent-by-stimulus matrix of dissimilarities.
#' @inheritParams smacf
#' @return A list of class `smunf`: `xr`, `xs`, `stress`, `iter` and
#'   `conv`.
#' @references Borg, I. and Groenen, P. J. F. (2005). *Modern
#'   Multidimensional Scaling*, 2nd ed. Springer.
#'   \doi{10.1007/0-387-28981-X}
#' @examples
#' set.seed(1)
#' smunf(matrix(runif(12), 3), maxit = 20)$stress
#' @family scaling
#' @export
smunf <- function(D, ndim = 2, maxit = 300, tol = 1e-6) {
  .sv("smunf", list(D = D, ndim = ndim, maxit = maxit, tol = tol))
}

#' Ideal points from an unfolding
#'
#' In an unfolding the respondents' coordinates are their ideal points
#' (Armstrong et al. 2021, section 4.5): the respondent coordinates as a
#' matrix.
#'
#' @param xr Respondent coordinates, as [mlsmu()] or [smunf()] give them.
#' @param xs The stimulus coordinates (not used).
#' @return A list of class `ipsrec` with the matrix `ideal`.
#' @references Armstrong, D. A., Bakker, R., Carroll, R., Hare, C.,
#'   Poole, K. T. and Rosenthal, H. (2021). *Analyzing Spatial Models of
#'   Choice and Judgment*, 2nd ed. Chapman & Hall/CRC.
#' @examples
#' set.seed(1)
#' u <- smunf(matrix(runif(12), 3), maxit = 20)
#' ipsrec(u$xr, u$xs)$ideal
#' @family scaling
#' @export
ipsrec <- function(xr, xs = NULL) .sv("ipsrec", list(xr = xr, xs = xs))

#' NOMINATE utilities and vote probabilities
#'
#' The Gaussian utility of Poole and Rosenthal (1985): legislator \eqn{i}
#' gains \eqn{U = \beta \exp(-\frac{1}{2} \sum_k w_k^2 (x_{ik} -
#' z_{jk})^2)} from the outcome \eqn{z_j}; the probability of a yea is
#' the logistic function of the yea utility minus the nay utility.
#'
#' @param x Legislators' ideal points (legislators by dimensions).
#' @param zy The yea outcomes (votes by dimensions).
#' @param zn The nay outcomes (votes by dimensions).
#' @param beta The signal-to-noise ratio.
#' @param w `NULL` for weights of 1, or one weight per dimension.
#' @return A list of class `nomutl`: `uyea`, `unay`, `udiff` and `prob`
#'   (legislators by votes).
#' @references Poole, K. T. and Rosenthal, H. (1985). A spatial model for
#'   legislative roll call analysis. *American Journal of Political
#'   Science*, 29(2), 357-384. \doi{10.2307/2111172}
#' @examples
#' set.seed(1)
#' nomutl(matrix(rnorm(8), 4), matrix(rnorm(6), 3), matrix(rnorm(6), 3))$prob
#' @family scaling
#' @export
nomutl <- function(x, zy, zn, beta = 15, w = NULL) {
  .sv("nomutl", list(x = x, zy = zy, zn = zn, beta = beta, w = w))
}

#' One NOMINATE vote probability
#'
#' The probability of a yea of [nomutl()] for one legislator and one
#' vote.
#'
#' @param x The ideal point (one value per dimension).
#' @param zy The yea outcome.
#' @param zn The nay outcome.
#' @inheritParams nomutl
#' @return A list of class `nompr` with `prob`.
#' @references Poole, K. T. and Rosenthal, H. (1985). A spatial model for
#'   legislative roll call analysis. *American Journal of Political
#'   Science*, 29(2), 357-384. \doi{10.2307/2111172}
#' @examples
#' nompr(0.1, 0.3, -0.3)$prob
#' @family scaling
#' @export
nompr <- function(x, zy, zn, beta = 15, w = NULL) {
  .sv("nompr", list(x = x, zy = zy, zn = zn, beta = beta, w = w))
}

#' NOMINATE log-likelihood and classification
#'
#' The log-likelihood of the observed votes under the probabilities of
#' [nomutl()] (clamped to \eqn{[10^{-10}, 1 - 10^{-10}]}) and `gmp`, the
#' share of votes the model classifies correctly.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @inheritParams nomutl
#' @return A list of class `nomll`: `loglik`, `gmp`, `correct` and
#'   `total`.
#' @references Poole, K. T. and Rosenthal, H. (1997). *Congress: A
#'   Political-Economic History of Roll Call Voting*. Oxford University
#'   Press.
#' @examples
#' set.seed(1)
#' v <- matrix(rbinom(20, 1, 0.5), 4, 5)
#' nomll(v, matrix(rnorm(4)), matrix(rnorm(5)), matrix(rnorm(5)))
#' @family scaling
#' @export
nomll <- function(votes, x, zy, zn, beta = 15, w = NULL) {
  .sv(
    "nomll",
    list(votes = votes, x = x, zy = zy, zn = zn, beta = beta, w = w)
  )
}

#' Procrustes rotation
#'
#' The orthogonal rotation (a proper rotation: no reflection) that best
#' carries the centred configuration `X` onto the centred `target` in
#' least squares, the target's centre added back (Gower and Dijksterhuis
#' 2004).
#'
#' @param X The configuration to rotate.
#' @param target The target configuration, of the same dimensions.
#' @return A list of class `procrs`: `rotated`, `rotation`, `scale` (the
#'   sum of the singular values of \eqn{X_c' T_c}) and `mse`.
#' @references Gower, J. C. and Dijksterhuis, G. B. (2004). *Procrustes
#'   Problems*. Oxford University Press.
#' @examples
#' set.seed(1)
#' A <- matrix(rnorm(20), 10)
#' procrs(A, A %*% matrix(c(0, 1, -1, 0), 2))$rotation
#' @family scaling
#' @export
procrs <- function(X, target) .sv("procrs", list(X = X, target = target))

#' Binary IRT likelihood
#'
#' The log-likelihood of the votes under the probit item response model
#' of Clinton, Jackman and Rivers (2004), \eqn{P(y_{ij} = 1) =
#' \Phi(\beta_j' x_i - \alpha_j)} (clamped to \eqn{[10^{-10}, 1 -
#' 10^{-10}]}), with the share of votes classified correctly.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @param x Ideal points (legislators by dimensions).
#' @param alpha One difficulty per vote.
#' @param beta Discriminations (votes by dimensions).
#' @return A list of class `irtlik`: `loglik`, `prob` (0 where absent),
#'   `correct`, `total` and `accuracy`.
#' @references Clinton, J., Jackman, S. and Rivers, D. (2004). The
#'   statistical analysis of roll call data. *American Political Science
#'   Review*, 98(2), 355-370. \doi{10.1017/S0003055404001194}
#' @examples
#' set.seed(1)
#' v <- matrix(rbinom(20, 1, 0.5), 4, 5)
#' irtlik(v, matrix(rnorm(4)), rep(0, 5), matrix(rnorm(5)))$loglik
#' @family scaling
#' @export
irtlik <- function(votes, x, alpha, beta) {
  .sv("irtlik", list(votes = votes, x = x, alpha = alpha, beta = beta))
}

#' Posterior summaries of ideal-point draws
#'
#' Means, standard deviations and central 95% intervals (R's default
#' quantile, type 7) of a chain of ideal points, each draw first
#' standardised dimension by dimension when `std` is `TRUE` (location and
#' scale are not identified in IRT models; Jackman 2009).
#'
#' @param chain An array of draws by legislators by dimensions (a matrix
#'   for one dimension).
#' @param std Whether to standardise each draw.
#' @return A list of class `irtpst`: `mean`, `sd`, `lo` and `hi`
#'   (legislators by dimensions), `nsamp` and `std`.
#' @references Jackman, S. (2009). *Bayesian Analysis for the Social
#'   Sciences*. Wiley.
#' @examples
#' set.seed(1)
#' irtpst(array(rnorm(100 * 5 * 2), c(100, 5, 2)))$mean
#' @family scaling
#' @export
irtpst <- function(chain, std = TRUE) {
  .sv("irtpst", list(chain = chain, std = std))
}

#' Ordered Optimal Classification
#'
#' Ordinal responses (a survey's agree-disagree items, a court's graded
#' outcomes) classified by item directions and cutpoints in the spirit of
#' Hare, Liu and Lupton (2018): respondents and item directions start
#' normal (seed 42, as [set.seed()] takes it), the cutpoints of each item
#' at equal steps on \eqn{[-1, 1]}; each sweep moves every cutpoint to
#' the median of the midpoints between neighbouring respondents (in
#' projection order) on either side of it, until the share classified
#' correctly stops changing.
#'
#' @param Y A respondent-by-item matrix of ordinal responses, `NA` for
#'   missing.
#' @param ndim The number of dimensions.
#' @param maxit The most sweeps.
#' @return A list of class `ordocl`: `ideal`, `cuts` (one vector per
#'   item), `normals`, `correct` (the share classified correctly) and
#'   `iter`.
#' @references Hare, C., Liu, T.-P. and Lupton, R. N. (2018). What
#'   Ordered Optimal Classification reveals about ideological structure,
#'   cleavages, and polarization in the American mass public. *Public
#'   Choice*, 176(1), 57-78.
#' @examples
#' set.seed(1)
#' Y <- matrix(sample(1:4, 60, replace = TRUE), 15, 4)
#' ordocl(Y, ndim = 1, maxit = 20)$correct
#' @family scaling
#' @export
ordocl <- function(Y, ndim = 2, maxit = 500) {
  .sv("ordocl", list(Y = Y, ndim = ndim, maxit = maxit))
}

#' Anchoring vignettes
#'
#' Corrects self-ratings for differential item functioning with
#' anchoring vignettes, nonparametrically (King et al. 2004): each
#' respondent's ratings of the vignettes, sorted, give that respondent's
#' thresholds between categories, and the self-rating is recoded as the
#' number of thresholds at or below it. The spread of each threshold
#' across respondents measures the differential item functioning.
#'
#' @param y Each respondent's self-rating, `NA` for missing.
#' @param V The respondents' ratings of the vignettes (respondents by
#'   vignettes), `NA` for missing.
#' @param ncat The number of ordered categories.
#' @return A list of class `anchv`: `scores` (the corrected ratings),
#'   `thresh` (respondents by thresholds), `dif` (the standard deviation
#'   of each threshold), `order` (the vignettes by mean rating), `n` and
#'   `nvig`.
#' @references King, G., Murray, C. J. L., Salomon, J. A. and Tandon, A.
#'   (2004). Enhancing the validity and cross-cultural comparability of
#'   measurement in survey research. *American Political Science Review*,
#'   98(1), 191-207. \doi{10.1017/S000305540400108X}
#' @examples
#' set.seed(1)
#' anchv(sample(1:5, 30, TRUE), matrix(sample(1:5, 90, TRUE), 30))$dif
#' @family scaling
#' @export
anchv <- function(y, V, ncat = 5) .sv("anchv", list(y = y, V = V, ncat = ncat))

#' INDSCAL: scaling with individual differences
#'
#' The weighted Euclidean model of Carroll and Chang (1970): one stimulus
#' space shared by every individual, each individual weighting its
#' dimensions. The configuration starts from classical scaling of the
#' mean dissimilarities; then the weights of each individual and
#' dimension (at least 0.01) and the stimulus coordinates are updated in
#' turn until the configuration changes by a relative `tol`. The stress
#' is the sum of squared differences between each individual's
#' dissimilarities and the distances in its weighted space.
#'
#' @param D A list of square dissimilarity matrices of one size, one per
#'   individual (judge, panel, survey wave).
#' @param ndim The number of dimensions.
#' @param maxit The most sweeps.
#' @param tol The convergence tolerance on the relative change of the
#'   configuration.
#' @return A list of class `indscl`: `config` (stimuli by dimensions),
#'   `weights` (individuals by dimensions), `stress`, `iter`, `nind` and
#'   `nstim`.
#' @references Carroll, J. D. and Chang, J.-J. (1970). Analysis of
#'   individual differences in multidimensional scaling via an N-way
#'   generalization of "Eckart-Young" decomposition. *Psychometrika*,
#'   35(3), 283-319. \doi{10.1007/BF02310791}
#' @examples
#' set.seed(1)
#' D1 <- as.matrix(dist(matrix(rnorm(20), 5)))
#' D2 <- as.matrix(dist(matrix(rnorm(20), 5)))
#' indscl(list(D1, D2), maxit = 30)$weights
#' @family scaling
#' @export
indscl <- function(D, ndim = 2, maxit = 300, tol = 1e-6) {
  .sv("indscl", list(D = D, ndim = ndim, maxit = maxit, tol = tol))
}

#' Normal vector of an outside measure
#'
#' Projects an outside measure (a judge's party, a legislator's vote
#' share) into an ideal-point space: the least-squares regression of the
#' measure on the coordinates, its slopes scaled to unit length give the
#' direction along which the measure rises, with its angle and the
#' \eqn{R^2} (Armstrong et al. 2021, section 2.6).
#'
#' @param x Ideal points (units by dimensions).
#' @param y The outside measure, `NA` for missing.
#' @return A list of class `nrmvec`: `normal`, `degrees`, `radians`,
#'   `r2` and `coef` (intercept first).
#' @references Armstrong, D. A., Bakker, R., Carroll, R., Hare, C.,
#'   Poole, K. T. and Rosenthal, H. (2021). *Analyzing Spatial Models of
#'   Choice and Judgment*, 2nd ed. Chapman & Hall/CRC.
#' @examples
#' set.seed(1)
#' nrmvec(matrix(rnorm(20), 10), rnorm(10))$degrees
#' @family scaling
#' @export
nrmvec <- function(x, y) .sv("nrmvec", list(x = x, y = y))

#' Cutting lines of votes in two dimensions
#'
#' The ends of each vote's cutting line \eqn{n_1 x + n_2 y = c} across
#' the horizontal range `xlim` (or from -10 to 10 vertically when the
#' line is vertical), its midpoint and the angle of its normal: the
#' lines of a Coombs mesh (Poole 2005). The plot itself is not drawn
#' here.
#'
#' @param normals The normal vectors of the votes (votes by 2).
#' @param cuts Each vote's cutpoint \eqn{c}.
#' @param xlim The horizontal range.
#' @return A list of class `cutlns`: `ends` (votes by 4: x1, y1, x2,
#'   y2), `angles` (degrees), `mids` (votes by 2) and `nlines`.
#' @references Poole, K. T. (2005). *Spatial Models of Parliamentary
#'   Voting*. Cambridge University Press.
#' @examples
#' set.seed(1)
#' cutlns(matrix(rnorm(6), 3), c(0.1, -0.2, 0))$ends
#' @family scaling
#' @export
cutlns <- function(normals, cuts, xlim = c(-1, 1)) {
  .sv("cutlns", list(normals = normals, cuts = cuts, xlim = xlim))
}

#' DW-NOMINATE-style ideal points
#'
#' A deterministic fixed-point fit in the Gaussian-utility family of
#' Poole and Rosenthal (1997): from normal starts (seed 42, as
#' [set.seed()] takes it), each vote's normal vector and midpoint follow
#' the mean positions of its two sides, and each legislator moves to the
#' mean of the outcome locations (midpoint plus or minus half the normal)
#' of the sides taken, the cloud shrunk into the unit ball, for `maxit`
#' sweeps. The log-likelihood and the share classified correctly use the
#' probit of the Gaussian utility difference with \eqn{\beta = 15} and
#' equal dimension weights.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @param ndim The number of dimensions.
#' @param maxit The sweeps.
#' @return A list of class `dwnom`: `ideal`, `weights`, `normals`, `mids`
#'   (the outcome midpoints), `cuts`, `loglik`, `gmp` and `ndim`.
#' @references Poole, K. T. and Rosenthal, H. (1997). *Congress: A
#'   Political-Economic History of Roll Call Voting*. Oxford University
#'   Press.
#' @examples
#' set.seed(1)
#' v <- matrix(rbinom(600, 1, 0.5), 20, 30)
#' dwnom(v, maxit = 20)$gmp
#' @family scaling
#' @export
dwnom <- function(votes, ndim = 2, maxit = 100) {
  .sv("dwnom", list(votes = votes, ndim = ndim, maxit = maxit))
}

#' Parametric bootstrap of ideal points
#'
#' Standard errors of ideal points by the parametric bootstrap of Lewis
#' and Poole (2004): votes are simulated from the fitted logistic
#' probabilities (\eqn{\beta = 15} times the projection minus the
#' cutpoint), absences kept, and each legislator re-estimated by 20
#' gradient steps from the fitted point plus normal noise of standard
#' deviation 0.1. Draws come from `seed`, as [set.seed()] takes it.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @param x Fitted ideal points.
#' @param normals Fitted normal vectors (votes by dimensions).
#' @param cuts Fitted cutpoints, one per vote.
#' @param nboot The bootstrap replicates.
#' @param seed The seed.
#' @return A list of class `nombt`: `se` and `mean` (legislators by
#'   dimensions) and `nboot`.
#' @references Lewis, J. B. and Poole, K. T. (2004). Measuring bias and
#'   uncertainty in ideal point estimates via the parametric bootstrap.
#'   *Political Analysis*, 12(2), 105-127.
#' @examples
#' set.seed(1)
#' v <- matrix(rbinom(40, 1, 0.5), 5, 8)
#' f <- dwnom(v, maxit = 5)
#' nombt(v, f$ideal, f$normals, f$cuts, nboot = 5)$se
#' @family scaling
#' @export
nombt <- function(votes, x, normals, cuts, nboot = 100, seed = 42) {
  .sv(
    "nombt",
    list(votes = votes, x = x, normals = normals, cuts = cuts,
         nboot = nboot, seed = seed)
  )
}

#' EM estimation of binary IRT ideal points
#'
#' The EM algorithm of Imai, Lo and Olmsted (2016) for the probit model
#' \eqn{y^*_{ij} = \alpha_j + \beta_j' x_i + e_{ij}}, yea when
#' \eqn{y^*_{ij} > 0}, with priors \eqn{x_i \sim N(\mu_x, s^2_x I)} and
#' \eqn{(\alpha_j, \beta_j) \sim N(\mu_b, s^2_b I)}. The E-step replaces
#' \eqn{y^*} by its truncated-normal mean (untruncated when absent); the
#' M-step is the pair of ridge regressions; the fixed point is the
#' posterior mode. The default start is \eqn{\alpha = \beta = 0} and
#' \eqn{x} from 50 steps of orthogonal iteration on the cross-product of
#' the column-centred +1/-1 votes.
#'
#' @param votes A legislator-by-vote matrix: 1 yea, 0 nay, `NA` absent.
#' @param ndim The number of dimensions.
#' @param maxit The most iterations.
#' @param tol The convergence threshold.
#' @param xprior The prior mean and variance of the ideal points.
#' @param bprior The prior mean and variance of the intercepts and
#'   slopes.
#' @param start `NULL`, or a list of the starting intercepts, slopes and
#'   ideal points.
#' @param conv `"cor"` (1 minus the smallest correlation of old and new
#'   values, column by column) or `"abs"` (the largest absolute change).
#' @return A list of class `irtem`: `ideal`, `disc`, `diff` (the
#'   intercepts), `loglik`, `iter` and `conv`.
#' @references Imai, K., Lo, J. and Olmsted, J. (2016). Fast estimation
#'   of ideal points with massive data. *American Political Science
#'   Review*, 110(4), 631-656. \doi{10.1017/S000305541600037X}
#' @examples
#' set.seed(1)
#' x <- rnorm(20)
#' v <- matrix(rbinom(600, 1, pnorm(outer(x, rnorm(30)))), 20)
#' cor(irtem(v)$ideal[, 1], x)
#' @family scaling
#' @export
irtem <- function(votes, ndim = 1, maxit = 500, tol = 1e-6,
                  xprior = c(0, 1), bprior = c(0, 25), start = NULL,
                  conv = "cor") {
  .sv(
    "irtem",
    list(votes = votes, ndim = ndim, maxit = maxit, tol = tol,
         xprior = xprior, bprior = bprior, start = start, conv = conv)
  )
}

#' Bootstrap of scaling positions
#'
#' Standard errors and percentile intervals of Aldrich-McKelvey stimulus
#' positions ([amscal()]) or Basic Space weights ([bbscal()], first
#' dimension, of the matrix or of its transpose) by resampling
#' respondents with replacement (Efron and Tibshirani 1993). Draws come
#' from `seed`, as [set.seed()] takes it. The sign of a Basic Space
#' dimension is arbitrary, so its replicates are best compared after
#' aligning them.
#'
#' @param Z A respondent-by-stimulus (or issue) matrix, `NA` for missing.
#' @param method `"am"`, `"blackbox"` or `"blackbox_t"`.
#' @param nboot The replicates.
#' @param seed The seed.
#' @return A list of class `scboot`: `se`, `mean`, `lo`, `hi` (2.5% and
#'   97.5%) and `nboot`.
#' @references Efron, B. and Tibshirani, R. J. (1993). *An Introduction to
#'   the Bootstrap*. Chapman & Hall.
#' @examples
#' set.seed(1)
#' Z <- matrix(rnorm(100), 20, 5) + rep(1:5, each = 20)
#' scboot(Z, nboot = 10)$se
#' @family scaling
#' @export
scboot <- function(Z, method = "am", nboot = 200, seed = 42) {
  .sv("scboot", list(Z = Z, method = method, nboot = nboot, seed = seed))
}

#' Wordfish: positions from word counts
#'
#' The Poisson scaling model of Slapin and Proksch (2008),
#' \eqn{y_{ij} \sim \mathrm{Poisson}(\exp(\psi_i + \alpha_j + \beta_j
#' \omega_i))}, for documents (judgments, briefs, platforms) by their
#' word counts: alternating Newton steps on the document positions
#' \eqn{\omega} (standardised after each sweep, the word parameters
#' adjusted to keep the linear predictor) and on the word parameters,
#' from normal starts (seed 42, as [set.seed()] takes it), until the
#' positions change by a relative `tol`. The sign of the scale is
#' arbitrary.
#'
#' @param dtm A document-by-term matrix of counts.
#' @param maxit The most sweeps.
#' @param tol The convergence tolerance.
#' @return A list of class `wrdfsh`: `pos` (\eqn{\omega}), `beta`,
#'   `alpha`, `psi`, `loglik` and `iter`.
#' @references Slapin, J. B. and Proksch, S.-O. (2008). A scaling model
#'   for estimating time-series party positions from texts. *American
#'   Journal of Political Science*, 52(3), 705-722.
#'   \doi{10.1111/j.1540-5907.2008.00338.x}
#' @examples
#' set.seed(1)
#' wrdfsh(matrix(rpois(600, 5), 20), maxit = 20)$pos
#' @family scaling
#' @export
wrdfsh <- function(dtm, maxit = 100, tol = 1e-6) {
  .sv("wrdfsh", list(dtm = dtm, maxit = maxit, tol = tol))
}
