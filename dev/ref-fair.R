# The rmorie values the fairness tests compare with, recorded in
# tests/testthat/ref-fair.rds. Run from the package root with rmorie
# installed (the owner's earlier package; its library LIB):
#   Rscript dev/ref-fair.R LIB
# The inputs are made here, from fixed seeds; the tests rebuild them with
# the same code (fairin() below, copied in test-fairdf.R).
a <- commandArgs(TRUE)
if (length(a)) .libPaths(c(a[1], .libPaths()))
suppressMessages(library(rmorie))
fairin <- function() {
  set.seed(2)
  n <- 300
  g <- sample(c("A", "B", "C"), n, TRUE)
  y <- rbinom(n, 1, 0.4)
  p <- rbinom(n, 1, c(A = 0.5, B = 0.35, C = 0.2)[g])
  per <- sample(c("p1", "p2", "p3"), n, TRUE)
  cty <- sample(c("X", "Y"), n, TRUE)
  list(g = g, y = y, p = p, per = per, cty = cty)
}
d <- fairin()
quiet <- function(e) suppressWarnings(suppressMessages(e))
r <- list(R = as.character(getRversion()),
          rmorie = as.character(packageVersion("rmorie")))
r$di <- quiet(morie_fairness_disparate_impact(d$p, d$g, privileged = "A"))$ratios
r$dp <- quiet(morie_fairness_demographic_parity(d$p, d$g, privileged = "A"))$gaps
b <- quiet(morie_fairness_bias_amplification(d$p, d$g, privileged = "A"))
r$bas <- c(bas = b$value, gap = b$demographic_parity_gap, gini = b$gini)
e <- quiet(morie_fairness_equalized_odds(d$y, d$p, d$g, privileged = "A"))
r$eo <- list(tpr = unlist(e$tpr_gaps), fpr = unlist(e$fpr_gaps),
             value = e$value)
ao <- quiet(morie_fairness_average_odds_difference(d$y, d$p, d$g,
                                                   privileged = "A"))
r$aod <- list(aod = unlist(ao$average_odds_difference), value = ao$value)
# no reference given: rmorie takes the highest share for the ratio, but the
# highest true positive rate for equalized odds
r$aodnull <- quiet(morie_fairness_average_odds_difference(d$y, d$p,
                                                          d$g))$privileged
r$di0 <- quiet(morie_fairness_disparate_impact(c(0, 0, 1, 0), c("a", "a",
                                                                "b", "b"),
                                               privileged = "a"))$ratios
tp <- quiet(morie_fairness_predpol_temporal_audit(d$per, d$cty, d$p, d$g,
                                                  privileged = "A"))$payload
r$tmp <- list(cells = do.call(rbind, lapply(tp$cells, function(x) {
  data.frame(x[c("city", "period")], n = x$n, dir = x$dir,
             pg = x$parity_gap, gini = x$gini, bas = x$bas)
})), city = tp$per_city, rng = tp$worst_dir_range,
spread = tp$cross_city_dir_spread)
tp1 <- quiet(morie_fairness_predpol_temporal_audit(d$per, rep("X", 300), d$p,
                                                   d$g, privileged = "A"))
r$spread1 <- tp1$payload$cross_city_dir_spread
r$bnd <- morie_fairness_base_rate_bounds(c(0.35, 0.55, 0.2), 0.1, 0.2)
r$cmp <- lapply(list(c(0.35, 0.55, 0.05, 0.20), c(0.35, 0.55, 0.05, 0.40),
                     c(0.6, 0.2, 0.1, 0.1), c(0.30, 0.31, 0, 0.05)),
                function(v) do.call(morie_fairness_compare_groups,
                                    as.list(v))[c("order",
                                                  "breakdown_beta_max")])
r$tr <- morie_fairness_true_rate(c(0.38, 0.5), 0.1, 0.2)
r$fpr <- morie_fairness_implied_fpr(0.3, 0.6, 0.25)
rk <- morie_ranking_resolution(c(0.2, 0.35, 0.8, 0.5), c(0.1, 0.1, 0.05, 0.2))
r$rank <- rk[c("identified_pairs", "share_unidentified", "resolution")]
v <- c(1, 2, 3, 4, 5, 6, 7, 8)
pm <- matrix(c(0, 1, 2, 4, 5, 7), ncol = 2)
r$lip <- Objfair(v, pm, L = 1.5, metric = c(1, 2, 0))[c(
  "L_required", "n_violations", "max_gap", "max_pair_i", "max_pair_j")]
ca <- quiet(morie_predpol_calibration_audit(
  areas = paste0("d", 1:8), mean_risk = c(90, 80, 70, 60, 30, 20, 10, 5),
  outcome_rate = c(10, 25, 20, 80, 70, 90, 85, 60),
  group = c("X", "X", "Y", "X", "Y", "Y", "X", "Y")))
r$cal <- ca[c("spearman", "spearman_pvalue", "group_rank_gap", "rank_gap")]
ca2 <- quiet(morie_fairness_predpol_calibration_audit(
  paste0("d", 1:8), c(90, 80, 70, 60, 30, 20, 10, 5),
  c(10, 25, 20, 80, 70, 90, 85, 60), c("X", "X", "Y", "X", "Y", "Y", "X",
                                       "Y")))
r$cal2 <- ca2[c("spearman", "spearman_pvalue")]
r$agg <- morie_predpol_aggregate_areas(
  area = c("a", "a", "b", "b", "c", "c", "c"),
  risk = c(10, 20, 30, 40, 5, 6, 7), outcome = c(1, 0, 1, 1, 0, 2, 1),
  group = c("x", "y", "y", "y", "x", "x", "y"),
  population = c(a = 100, b = 200, c = 50))
x7 <- c(9, 10, 11, 19, 20, 21, 30, 14)
g7 <- c("A", "A", "A", "B", "B", "B", "C", "A")
r$sd <- quiet(morie_predpol_score_disparity(x7, g7))[c(
  "spread", "group_means", "gaps", "anova_f", "anova_pvalue", "reference")]
r$sd2 <- quiet(morie_fairness_predpol_score_disparity(x7, g7))[c(
  "anova_f", "anova_pvalue")]
set.seed(1)
cx <- matrix(runif(20), ncol = 2L)
ox <- matrix(runif(10), ncol = 2L)
nd <- morie_fairness_noisy_or_detection(cx, ox, radius = 0.5, p_detect = 0.7)
r$nor <- list(cx = cx, ox = ox, prob = nd$probabilities,
              k = nd$officers_in_range)
set.seed(3)
n <- 400
g <- rep(c("A", "B"), each = n / 2)
sc <- rnorm(n, ifelse(g == "A", 5, 4))
pr <- rpois(n, ifelse(g == "A", 2, 1))
sen <- 2 + 0.8 * sc + 0.5 * pr + ifelse(g == "A", 1, 0) + rnorm(n)
od <- data.frame(sen, sc, pr, g)
o <- morie_disparity_decomposition(sen ~ sc + pr, od, group = "g",
                                   reference = "B")
r$oax <- list(d = od, o = o[c("gap", "twofold_B", "twofold_A", "threefold",
                              "by_variable")])
saveRDS(r, "tests/testthat/ref-fair.rds", version = 2)
cat("recorded with R", r$R, "and rmorie", r$rmorie, "\n")
