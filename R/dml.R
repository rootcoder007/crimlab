#' Double machine learning by cross-fitting
#'
#' Estimates a causal parameter from orthogonal scores whose nuisance
#' functions (outcome regressions, propensities) are learned out of fold:
#' each fold's predictions come from learners fitted on the other folds,
#' so that the estimate keeps root-n inference while the learners are
#' flexible (Chernozhukov et al. 2018). The learning, the scores, the
#' estimate and its standard error are all compiled; the folds are the
#' caller's, so results can be matched with other software fold for fold.
#'
#' The models are those of 'DoubleML', with its scores and its defaults:
#' * `"plr"`, the partially linear regression `y = theta d + g(x) + e`,
#'   score `"po"` (partialling out: nuisances `l = E[y | x]`, `m = E[d |
#'   x]`) or `"iv"` (IV-type: also `g`, learned on `y - theta0 d` with
#'   `theta0` from the partialling-out score);
#' * `"pliv"`, the partially linear IV model with one instrument `z`
#'   (nuisances `l`, `m = E[z | x]`, `r = E[d | x]`, and `g` for `"iv"`);
#' * `"irm"`, the interactive regression model for a 0/1 treatment, score
#'   `"ate"` or `"atte"` (nuisances `g0`, `g1`, the outcome regressions
#'   among the untreated and the treated, and the propensity `m`);
#' * `"iivm"`, the interactive IV model (local average effect) for a 0/1
#'   treatment and a 0/1 instrument (`g0`, `g1`, `m = P(z = 1 | x)`, `r0`,
#'   `r1`, the treatment probabilities given `z`);
#' * `"did"`, difference in differences on a panel: `y` is the change in
#'   the outcome, `d` the 0/1 group (nuisances `g0`, the change among the
#'   untreated, and `m`), with the in-sample normalised observational score
#'   of 'DoubleML' for 'Python' (Sant'Anna and Zhao 2020; Chang 2020);
#' * `"didcs"`, the same on repeated cross-sections, `t` the 0/1 period
#'   (nuisances `g00`, `g01`, `g10`, `g11`, the outcome regressions in each
#'   group and period, and `m`);
#' * `"ssm"`, the average effect when the outcome is observed only where
#'   `s` is 1, missing at random given `d` and `x` (Bia, Huber and Laffers
#'   2023; nuisances `pi = P(s = 1 | d, x)`, `m`, `g0`, `g1`);
#' * `"aipw"` and `"tmle"`, the average effect of a 0/1 treatment from the
#'   nuisances of `"irm"`: augmented inverse probability weighting, with the
#'   standard error of the 'AIPW' package (the standard deviation of the
#'   influence function), and targeted maximum likelihood with the
#'   fluctuation of the 'tmle' package (Gruber and van der Laan 2012).
#'
#' The estimate solves the mean of the score over all folds (DoubleML's
#' `"dml2"`); its variance is `mean(psi^2) / mean(psi_a)^2 / n`, or, with
#' clusters (one-way, models `"plr"` and `"pliv"`), the cluster-robust form
#' of Chiang et al. (2022). With several repetitions of the folds (columns
#' of `fold`), the estimate is the median over repetitions and the
#' variance the median of `var + (theta_r - theta)^2`.
#'
#' @param y The outcome: one number per row of `x` (for `"did"`, the
#'   change between the periods; for `"ssm"`, `NA` allowed where `s` is 0).
#' @param d The treatment: numbers (`"plr"`, `"pliv"`) or 0 and 1.
#' @param x The covariates: a numeric matrix, one row per observation.
#' @param model One of `"plr"`, `"pliv"`, `"irm"`, `"iivm"`, `"aipw"`,
#'   `"tmle"`, `"did"`, `"didcs"`, `"ssm"`.
#' @param fold The fold of each row, whole numbers from 1 to the number of
#'   folds (at least 2), each used; or a matrix of them, one column per
#'   repetition of the cross-fitting.
#' @param z The instrument (`"pliv"`; 0 and 1 for `"iivm"`), else `NULL`.
#' @param s For `"ssm"`, 1 where the outcome is observed and 0 where not.
#' @param t For `"didcs"`, the period of each row, 0 or 1.
#' @param score `NULL` for the default (`"po"`, or `"ate"`), or `"po"` or
#'   `"iv"` for `"plr"` and `"pliv"`, `"ate"` or `"atte"` for `"irm"`.
#' @param lrn The learner of the outcome nuisances: `"glm"` (least squares;
#'   logistic regression for a probability), `"lasso"`, `"ridge"` or
#'   `"enet"` (the elastic net of [enpath()] with mixing `alpha`, the penalty
#'   chosen by cross-validation within the training folds), or `"forest"`
#'   (the random forest of [rfor()] with `ntree` trees and `mtry` columns
#'   per split, nodes of 5 rows, or 10 for a probability, left unsplit).
#' @param lrnd The learner of the treatment, instrument and selection
#'   nuisances, from the same choices.
#' @param alpha The mixing parameter of `"enet"`, from 0 to 1.
#' @param nlam The number of penalties on each elastic net path.
#' @param nfin The number of inner folds that choose the penalty, at
#'   least 2.
#' @param ntree The number of trees of `"forest"`: a whole number, at least 1.
#' @param mtry The columns tried at each split of `"forest"`: `NULL` for the
#'   floor of the square root of the number of columns, or a whole number
#'   (more than the columns means all of them).
#' @param trim Propensities are kept within `trim` of 0 and 1: `NULL` for
#'   the default of the reference implementation (1e-12 for `"irm"`,
#'   `"iivm"` and `"ssm"`, as 'DoubleML'; 0.01 for `"did"` and `"didcs"`, as
#'   'DoubleML' for 'Python'; 0.025 for `"aipw"`, as 'AIPW'; for `"tmle"`,
#'   `5 / sqrt(n) / log(n)`, a lower bound only, as 'tmle').
#' @param clu `NULL`, or the cluster of each row (whole numbers from 1 to
#'   the number of clusters, each used), for `"plr"` and `"pliv"`; folds
#'   must then be constant within clusters.
#' @param pred `NULL` to learn the nuisances, or their predictions: a
#'   matrix with one row per observation and one column per nuisance, in
#'   the order listed above (an array with a third dimension for several
#'   repetitions). The estimate is then a closed form of these.
#' @param seed The seed of the inner folds and of the forests: one number,
#'   at least 0.
#' @return A list: `coef` and `se`, the estimate and its standard error;
#'   `reps`, a data frame with each repetition's `coef` and `se`; `psia`
#'   and `psib`, the score elements (one column per repetition; the
#'   estimate solves `mean(psia * coef + psib) = 0`); `pred`, the
#'   out-of-fold predictions, an array named by nuisance.
#' @references Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E.,
#'   Hansen, C., Newey, W. and Robins, J. (2018). Double/debiased machine
#'   learning for treatment and structural parameters. *The Econometrics
#'   Journal*, 21(1), C1-C68. \doi{10.1111/ectj.12097}
#'
#'   Chiang, H. D., Kato, K., Ma, Y. and Sasaki, Y. (2022). Multiway
#'   cluster robust double/debiased machine learning. *Journal of Business
#'   and Economic Statistics*, 40(3), 1046-1056.
#'   \doi{10.1080/07350015.2021.1895815}
#'
#'   Sant'Anna, P. H. C. and Zhao, J. (2020). Doubly robust
#'   difference-in-differences estimators. *Journal of Econometrics*,
#'   219(1), 101-122. \doi{10.1016/j.jeconom.2020.06.003}
#'
#'   Chang, N.-C. (2020). Double/debiased machine learning for
#'   difference-in-differences models. *The Econometrics Journal*, 23(2),
#'   177-191. \doi{10.1093/ectj/utaa001}
#'
#'   Bia, M., Huber, M. and Laffers, L. (2023). Double machine learning for
#'   sample selection models. *Journal of Business and Economic
#'   Statistics*, 42(3), 958-969. \doi{10.1080/07350015.2023.2271071}
#'
#'   Gruber, S. and van der Laan, M. J. (2012). tmle: an R package for
#'   targeted maximum likelihood estimation. *Journal of Statistical
#'   Software*, 51(13), 1-35. \doi{10.18637/jss.v051.i13}
#' @examples
#' set.seed(1)
#' n <- 500
#' x <- matrix(rnorm(n * 5), n)
#' d <- x[, 1] + rnorm(n)
#' y <- 0.5 * d + x[, 1] - x[, 2] + rnorm(n)
#' fold <- sample(rep(1:5, length.out = n))
#' r <- dmlx(y, d, x, "plr", fold, lrn = "glm", lrnd = "glm")
#' c(r$coef, r$se)
#' @export
dmlx <- function(y, d, x, model = "plr", fold, z = NULL, s = NULL, t = NULL,
                 score = NULL, lrn = "lasso", lrnd = lrn, alpha = 0.5,
                 nlam = 100, nfin = 5, ntree = 500, mtry = NULL, trim = NULL,
                 clu = NULL, pred = NULL, seed = 1) {
  .Call(C_clab_r_dmlx, y, d, x, model, fold, z, s, t, score, lrn, lrnd,
        alpha, nlam, nfin, ntree, mtry, trim, clu, pred, seed)
}
