# Bayesian scaling and ideal points: thin wrappers over the samplers of
# src/svb.cpp (src/svk.h).

#' Bayesian Aldrich-McKelvey scaling
#'
#' Respondent \eqn{i} places stimulus \eqn{j} at \eqn{z_{ij} \sim N(a_i +
#' b_i \hat z_j, 1 / (\tau_i \tau_j))}: a respondent's shift and stretch
#' of the common positions \eqn{\hat z}, with noise whose precision
#' factors into a respondent and a stimulus part (Hare et al. 2015).
#' Priors: \eqn{a_i, b_i \sim U(-100, 100)}, \eqn{\tau_j \sim G(0.1,
#' 0.1)}, \eqn{\tau_i \sim G(g_a, g_b)} with \eqn{g_a, g_b \sim G(0.1,
#' 0.1)}, and \eqn{\hat z} the standardised \eqn{z^* \sim N(0, 1)}, the
#' `polarity` stimulus restricted to the left. Gibbs sampling:
#' \eqn{a_i} and \eqn{b_i} from truncated normals, the precisions and
#' \eqn{g_b} from conjugate gammas, \eqn{g_a} and each \eqn{z^*_j} by
#' slice sampling (Neal 2003). Draws come from `seed`, as [set.seed()]
#' takes it; rows with no placement are dropped.
#'
#' @param Z A respondent-by-stimulus matrix of placements, `NA` for
#'   missing.
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` must be at
#'   most 2147483647.
#' @param polarity The column of the stimulus placed on the left.
#' @param seed The seed.
#' @return A list of class `bamscl`: `zeta`, `sd` and `interval` (2.5%
#'   and 97.5% rows) of the stimulus positions, `a` and `b` (posterior
#'   means), `tau` (stimulus precisions), `draws` (draws by stimuli) and
#'   `nsamp`.
#' @references Hare, C., Armstrong, D. A., Bakker, R., Carroll, R. and
#'   Poole, K. T. (2015). Using Bayesian Aldrich-McKelvey scaling to study
#'   citizens' ideological preferences and perceptions. *American Journal
#'   of Political Science*, 59(3), 759-774. \doi{10.1111/ajps.12151}
#'
#'   Neal, R. M. (2003). Slice sampling. *Annals of Statistics*, 31(3),
#'   705-767. \doi{10.1214/aos/1056562461}
#' @examples
#' # 40 simulated respondents place four stimuli
#' set.seed(1)
#' zt <- c(-1, -0.3, 0.4, 1.2)
#' Z <- rnorm(40, 0, 0.5) + outer(runif(40, 0.5, 1.5), zt) +
#'   matrix(rnorm(160, 0, 0.3), 40, 4)
#' round(bamscl(Z, nsamp = 100, burn = 50)$zeta, 2)
#' @family scaling
#' @export
bamscl <- function(Z, nsamp = 1000, burn = 200, polarity = 1, seed = 42) {
  .sv(
    "bamscl",
    list(Z = Z, nsamp = nsamp, burn = burn, polarity = polarity, seed = seed)
  )
}

#' Bayesian metric multidimensional scaling
#'
#' Each observed dissimilarity is lognormal about the distance between
#' the two objects' coordinates, \eqn{\log \delta_{ij} \sim N(\log
#' d_{ij}, 1 / \tau)}, with \eqn{x_{ik} \sim N(0, 10^2)} and \eqn{\tau
#' \sim U(0, 10)} (Bakker and Poole 2013). Coordinates are slice sampled
#' one at a time (Neal 2003), \eqn{\tau} drawn from its truncated gamma;
#' the chain starts from classical scaling plus small normal noise. The
#' likelihood is unchanged by translation, rotation and reflection, so
#' every draw is aligned (translation and rotation) onto the posterior
#' mean. Zero, negative or missing dissimilarities are left out. Draws
#' come from `seed`, as [set.seed()] takes it.
#'
#' @param D A square matrix of dissimilarities.
#' @inheritParams optcls
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` must be at
#'   most 2147483647.
#' @param sigma The starting \eqn{\sigma = \tau^{-1/2}}.
#' @param seed The seed.
#' @return A list of class `bmdsf`: `pos` (aligned posterior means), `sd`,
#'   `dist` (posterior mean distances), `sigma`, `tau`, `draws` (draws by
#'   objects by dimensions) and `nsamp`.
#' @references Bakker, R. and Poole, K. T. (2013). Bayesian metric
#'   multidimensional scaling. *Political Analysis*, 21(1), 125-140.
#'   \doi{10.1093/pan/mps039}
#' @examples
#' set.seed(1)
#' D <- as.matrix(dist(matrix(rnorm(16), 8)))
#' f <- bmdsf(D, nsamp = 50, burn = 50)
#' cor(f$dist[upper.tri(D)], D[upper.tri(D)])
#' @family scaling
#' @export
bmdsf <- function(D, ndim = 2, nsamp = 1000, burn = 200, sigma = 1,
                  seed = 42) {
  .sv(
    "bmdsf",
    list(D = D, ndim = ndim, nsamp = nsamp, burn = burn, sigma = sigma,
         seed = seed)
  )
}

#' Bayesian unfolding
#'
#' The lognormal model of [bmdsf()] for a respondent-by-stimulus matrix:
#' \eqn{\log \delta_{ij} \sim N(\log \lVert x_i - z_j \rVert, 1 / \tau)},
#' \eqn{N(0, 10^2)} priors on every coordinate and \eqn{\tau \sim U(0,
#' 10)} (Bakker and Poole 2013). Feeling thermometers \eqn{T} on 0-100
#' enter as \eqn{\delta = (100 - T) / 50}. Given the stimuli the
#' respondents are independent, and the other way round, so each block is
#' slice sampled at once, dimension by dimension; draws are aligned on
#' the stimuli. The chain starts from the singular value decomposition of
#' the centred matrix. Draws come from `seed`, as [set.seed()] takes it.
#'
#' @inheritParams mlsmu
#' @inheritParams bmdsf
#' @return A list of class `bunfld`: `stim` (aligned posterior means),
#'   `sd`, `ideal`, `dist` (posterior mean respondent-stimulus distances),
#'   `sigma`, `tau` and `nsamp`.
#' @references Bakker, R. and Poole, K. T. (2013). Bayesian metric
#'   multidimensional scaling. *Political Analysis*, 21(1), 125-140.
#'   \doi{10.1093/pan/mps039}
#' @examples
#' set.seed(1)
#' S <- matrix(runif(8, -0.5, 0.5), 4)
#' R <- matrix(runif(60, -0.5, 0.5), 30)
#' D <- sqrt(outer(R[, 1], S[, 1], "-")^2 + outer(R[, 2], S[, 2], "-")^2)
#' f <- bunfld(D, nsamp = 50, burn = 50)
#' cor(as.numeric(f$dist), as.numeric(D))
#' @family scaling
#' @export
bunfld <- function(D, ndim = 2, nsamp = 1000, burn = 200, seed = 42) {
  .sv(
    "bunfld",
    list(D = D, ndim = ndim, nsamp = nsamp, burn = burn, seed = seed)
  )
}

#' Clinton-Jackman-Rivers Bayesian IRT
#'
#' The probit ideal-point model \eqn{P(y_{ij} = 1) = \Phi(\beta_j' x_i -
#' \alpha_j)} with priors \eqn{x_i \sim N(0, I)} and \eqn{(\alpha_j,
#' \beta_j) \sim N(0, v I)} (Clinton, Jackman and Rivers 2004), by the
#' Gibbs sampler of Albert and Chib (1993): latent utilities from
#' truncated normals, then every \eqn{(\alpha_j, \beta_j)} and every
#' \eqn{x_i} from its normal full conditional. Every draw is an inverse
#' normal of a uniform of Philox4x32-10 (Salmon et al. 2011), sweep
#' \eqn{t} on stream \eqn{t} of key `seed`, so a chain is the same on
#' every platform. The default start is the posterior mode of [irtem()].
#'
#' @inheritParams optcls
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` must be at
#'   most 2147483647.
#' @param bvar The prior variance \eqn{v} of the vote parameters.
#' @param start `NULL`, or starting ideal points.
#' @param seed The Philox key, a whole number from 0.
#' @param chains Whether to return the draws as well.
#' @return A list of class `cjrirt`: `ideal`, `sd`, `alpha`, `beta`,
#'   `nsamp` and, with `chains`, `xchain` (draws by legislators by
#'   dimensions), `achain` and `bchain`.
#' @references Clinton, J., Jackman, S. and Rivers, D. (2004). The
#'   statistical analysis of roll call data. *American Political Science
#'   Review*, 98(2), 355-370. \doi{10.1017/S0003055404001194}
#'
#'   Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
#'   polychotomous response data. *Journal of the American Statistical
#'   Association*, 88(422), 669-679. \doi{10.1080/01621459.1993.10476321}
#'
#'   Salmon, J. K., Moraes, M. A., Dror, R. O. and Shaw, D. E. (2011).
#'   Parallel random numbers: as easy as 1, 2, 3. *Proceedings of the
#'   International Conference for High Performance Computing, Networking,
#'   Storage and Analysis* (SC11). \doi{10.1145/2063384.2063405}
#' @examples
#' votes <- rbind(c(1, 1, 0, 0), c(1, 1, 1, 0), c(0, 1, 1, 1),
#'                c(0, 0, 1, 1), c(1, 0, 0, 0))
#' round(cjrirt(votes, nsamp = 50, burn = 10)$ideal, 3)
#' @family scaling
#' @export
cjrirt <- function(votes, ndim = 1, nsamp = 1000, burn = 200, bvar = 25,
                   start = NULL, seed = 0, chains = TRUE) {
  .sv(
    "cjrirt",
    list(votes = votes, ndim = ndim, nsamp = nsamp, burn = burn,
         bvar = bvar, start = start, seed = seed, chains = chains)
  )
}

#' Alpha-NOMINATE
#'
#' Legislator \eqn{i} votes yea on roll call \eqn{j} with probability
#' \eqn{\Phi(u_{ij})}, \eqn{u_{ij} = Q_{ij} + \alpha (G_{ij} - Q_{ij})},
#' where \eqn{Q} is the quadratic and \eqn{G} the Gaussian utility
#' difference of the distances to the yea and nay outcomes, with weight
#' 0.5 (Carroll et al. 2013): \eqn{\alpha = 1} is Gaussian utility,
#' \eqn{\alpha = 0} quadratic. \eqn{\beta} has a flat prior on
#' \eqn{(0, \infty)}, \eqn{\alpha} on \eqn{[0, 1]}; coordinates carry the
#' hierarchical penalty \eqn{-x'Sx/2} with \eqn{S} redrawn each sweep
#' from an inverse Wishart. Every coordinate, \eqn{\beta} and
#' \eqn{\alpha} are slice sampled (Neal 2003: width 8, at most 3 steps
#' out); starts are uniform on \eqn{(-1, 1)}, \eqn{\beta = 10},
#' \eqn{\alpha = 0.7}. Each draw is then centred on its legislators and
#' rotated onto the posterior mean, the `polarity` legislator positive.
#' Roll calls whose minority side is below `lop` and legislators with
#' fewer than `minvotes` votes are dropped first. Draws come from `seed`,
#' as [set.seed()] takes it.
#'
#' @inheritParams optcls
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` times
#'   `thin` must be at most 2147483647.
#' @param seed The seed.
#' @param thin Keep every `thin`-th sweep after `burn`.
#' @param lop The smallest minority share of a kept roll call.
#' @param minvotes The fewest kept roll calls of a kept legislator.
#' @param polarity The row of the legislator placed on the positive side
#'   of each dimension (recycled).
#' @param fixalpha Whether to fix \eqn{\alpha = 1}.
#' @return A list of class `anomin`: `ideal` and `sd` (kept legislators
#'   by dimensions), `yea` and `nay` (outcome locations), `alpha`,
#'   `interval` (2.5% and 97.5% of \eqn{\alpha}), `beta`, `draws` (`X`,
#'   `Y`, `N` arrays, draws first, and `beta`, `alpha`), `legis` and
#'   `votes` (the kept rows and columns) and `ndim`.
#' @references Carroll, R., Lewis, J. B., Lo, J., Poole, K. T. and
#'   Rosenthal, H. (2013). The structure of utility in spatial models of
#'   voting. *American Journal of Political Science*, 57(4), 1008-1028.
#'   \doi{10.1111/ajps.12029}
#'
#'   Neal, R. M. (2003). Slice sampling. *Annals of Statistics*, 31(3),
#'   705-767. \doi{10.1214/aos/1056562461}
#' @examples
#' # 20 simulated legislators on a line, 30 roll calls, some votes flipped
#' set.seed(1)
#' x <- seq(-1, 1, length.out = 20)
#' V <- outer(x, runif(30, -0.7, 0.7), ">") * 1
#' flip <- sample(length(V), 30)
#' V[flip] <- 1 - V[flip]
#' f <- anomin(V, nsamp = 20, burn = 10, minvotes = 10, polarity = 20)
#' f$alpha
#' @family scaling
#' @export
anomin <- function(votes, ndim = 1, nsamp = 500, burn = 100, seed = 42,
                   thin = 1, lop = 0.025, minvotes = 20, polarity = 1,
                   fixalpha = FALSE) {
  .sv(
    "anomin",
    list(votes = votes, ndim = ndim, nsamp = nsamp, burn = burn,
         seed = seed, thin = thin, lop = lop, minvotes = minvotes,
         polarity = polarity, fixalpha = fixalpha)
  )
}

#' Ordinal IRT (ordinal factor model)
#'
#' Each ordinal response is a thresholded latent normal, \eqn{y^*_{ij} =
#' \lambda_{j0} + \lambda_j' \phi_i + e_{ij}}, \eqn{y_{ij} = c} when
#' \eqn{\gamma_{j,c-1} < y^*_{ij} \le \gamma_{j,c}}, with
#' \eqn{\gamma_{j1} = 0}, \eqn{\phi_i \sim N(0, I)}, loadings
#' \eqn{\lambda_j \sim N(0, I / L_0)} and a flat prior on the free
#' cutpoints (Quinn 2004). Each sweep updates the free cutpoints of every
#' item by the Metropolis-Hastings step of Cowles (1996) with \eqn{y^*}
#' integrated out (step sizes tuned during the burn-in only), then draws
#' \eqn{y^*}, the loadings and the ideal points from their full
#' conditionals. The chain starts from the leading principal components
#' of the mean-filled responses. Each draw is reflected so the first item
#' loads positively (rotated onto the first draw with several
#' dimensions). Draws come from `seed`, as [set.seed()] takes it.
#'
#' @param Y A respondent-by-item matrix of ordinal responses (any ordered
#'   numeric coding), `NA` for missing.
#' @inheritParams optcls
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` must be at
#'   most 2147483647.
#' @param seed The seed.
#' @param lprec The prior precision \eqn{L_0} of the loadings and
#'   intercepts (0, flat).
#' @return A list of class `ordirt`: `ideal` and `sd`, `disc` (items by
#'   dimensions), `icept`, `cuts` (one vector per item, the first 0),
#'   `accept` (cutpoint acceptance after the burn-in) and `nsamp`.
#' @references Quinn, K. M. (2004). Bayesian factor analysis for mixed
#'   ordinal and continuous responses. *Political Analysis*, 12(4),
#'   338-353.
#'
#'   Cowles, M. K. (1996). Accelerating Monte Carlo Markov chain
#'   convergence for cumulative-link generalized linear models.
#'   *Statistics and Computing*, 6(2), 101-111.
#' @examples
#' set.seed(1)
#' th <- rnorm(60)
#' Y <- matrix(cut(outer(th, runif(6, 0.8, 1.5)) + rnorm(360),
#'                 c(-Inf, -0.5, 0.5, Inf), labels = FALSE), 60, 6)
#' f <- ordirt(Y, nsamp = 50, burn = 50)
#' cor(f$ideal[, 1], th)
#' @family scaling
#' @export
ordirt <- function(Y, ndim = 1, nsamp = 500, burn = 100, seed = 42,
                   lprec = 0) {
  .sv(
    "ordirt",
    list(Y = Y, ndim = ndim, nsamp = nsamp, burn = burn, seed = seed,
         lprec = lprec)
  )
}

#' Dynamic IRT (ideal points that drift)
#'
#' The one-dimensional dynamic item response model of Martin and Quinn
#' (2002), made for courts: a vote is cast when \eqn{z_{jk} = -\alpha_k +
#' \beta_k \theta_{j,t(k)} + \varepsilon_{jk} > 0}; the ideal points
#' follow a random walk, \eqn{\theta_{j,0} \sim N(e_0, E_0)},
#' \eqn{\theta_{j,t} \sim N(\theta_{j,t-1}, \tau^2_j)}; \eqn{\alpha_k
#' \sim N(a_0, 1/A_0)}, \eqn{\beta_k \sim N(b_0, 1/B_0)}. Each sweep
#' draws the latent utilities, every case's \eqn{(\alpha, \beta)}, every
#' justice's path by forward filtering and backward sampling and, when
#' `c0` and `d0` are positive, the evolution variances from
#' \eqn{IG((c_0 + T)/2, (d_0 + SS)/2)}. A draw is reflected whenever the
#' `anchor` member's mean ideal point is negative. The chain starts from
#' the first principal component of the votes. Draws come from `seed`, as
#' [set.seed()] takes it.
#'
#' @param votes A member-by-vote matrix of 0/1 votes (-1/1 is read as
#'   0/1), `NA` absent.
#' @param period The period (term) of each vote, one per column.
#' @param nsamp The draws kept; `nsamp` times the values kept per draw
#'   must be at most 67108864 (2^26).
#' @param burn The sweeps discarded first; `burn` plus `nsamp` times
#'   `thin` must be at most 2147483647.
#' @param seed The seed.
#' @param thin Keep every `thin`-th sweep.
#' @param tau2 The evolution variance (its start when estimated).
#' @param c0,d0 The inverse-gamma prior of \eqn{\tau^2}, held at `tau2`
#'   unless both are positive.
#' @param e0,E0 The prior mean and variance of the first ideal point.
#' @param a0,A0 The prior mean and precision of \eqn{\alpha}.
#' @param b0,B0 The prior mean and precision of \eqn{\beta}.
#' @param anchor `NULL` (the largest first principal component score) or
#'   the row kept on the positive side.
#' @return A list of class `dynirt`: `theta` and `sd` (members by
#'   periods), `alpha`, `beta`, `tau2`, `periods`, `nsamp` and `anchor`.
#' @references Martin, A. D. and Quinn, K. M. (2002). Dynamic ideal point
#'   estimation via Markov chain Monte Carlo for the U.S. Supreme Court,
#'   1953-1999. *Political Analysis*, 10(2), 134-153.
#'   \doi{10.1093/pan/10.2.134}
#' @examples
#' # 12 simulated justices over three terms of 15 cases
#' set.seed(7)
#' period <- rep(1:3, each = 15)
#' walk <- function(s) s + cumsum(c(0, rnorm(2, 0, 0.3)))
#' truth <- t(sapply(rnorm(12), walk))
#' p <- pnorm(-rep(rnorm(45, 0, 0.5), each = 12) +
#'   rep(rnorm(45, 1.5, 0.5), each = 12) * truth[, period])
#' V <- matrix(rbinom(length(p), 1, p), 12)
#' f <- dynirt(V, period, nsamp = 100, burn = 50)
#' cor(as.numeric(f$theta), as.numeric(truth))
#' @family scaling
#' @export
dynirt <- function(votes, period, nsamp = 500, burn = 100, seed = 42,
                   thin = 1, tau2 = 1, c0 = -1, d0 = -1, e0 = 0, E0 = 1,
                   a0 = 0, A0 = 0.1, b0 = 0, B0 = 0.1, anchor = NULL) {
  .sv(
    "dynirt",
    list(votes = votes, period = period, nsamp = nsamp, burn = burn,
         seed = seed, thin = thin, tau2 = tau2, c0 = c0, d0 = d0, e0 = e0,
         E0 = E0, a0 = a0, A0 = A0, b0 = b0, B0 = B0, anchor = anchor)
  )
}
