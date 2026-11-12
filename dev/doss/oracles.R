# Reference values from the "oracle" packages for crimlab's existing
# forensic kernels (clab_gtmp, clab_kinidx, clab_patidx, clab_mixp,
# clab_cllr, clab_lrcal), for src/t_orc.cpp. R 4.5.0, with the user's
# ~/.Renviron (DNAtools in the dossier library needs Rcpp >= 1.1.0, which
# only the user library on zeus has).
# DNAtools first: it needs a newer Rcpp than the other packages load
source("dev/doss/pre.R")
suppressMessages({library(DNAtools); library(Familias); library(numberofalleles); library(comparison)})
for (p in c("Familias", "numberofalleles", "comparison")) cat(p, format(packageVersion(p)), "\n")

# one locus, alleles A..F with these frequencies, no mutation
fr <- c(0.05, 0.1, 0.15, 0.2, 0.22, 0.28); al <- LETTERS[1:6]
loc <- FamiliasLocus(frequencies = fr, allelenames = al, name = "L")
lr <- function(peds, data, kin = 0) {
  r <- FamiliasPosterior(peds, list(loc), data, ref = length(peds), kinship = kin)
  as.numeric(r$likelihoods)
}

# 1. paternity index: mother, child, alleged father; father vs unrelated
persons <- c("mo", "ch", "af")
pf <- FamiliasPedigree(id = persons, dadid = c(NA, "af", NA), momid = c(NA, "mo", NA),
                       sex = c("female", "female", "male"))
pu <- FamiliasPedigree(id = c(persons, "tf"), dadid = c(NA, "tf", NA, NA), momid = c(NA, "mo", NA, NA),
                       sex = c("female", "female", "male", "male"))
trios <- list(c("A","B", "A","C", "C","D"), c("A","B", "A","B", "B","B"), c("A","B", "A","B", "A","B"),
              c("C","C", "C","D", "D","E"), c("A","B", "B","C", "C","C"), c("E","F", "E","F", "A","B"),
              c("D","D", "D","D", "D","F"), c("B","C", "B","D", "D","D"))
cat("\n// Familias: LR(af father of ch | mo mother) / LR(unrelated), trios (mo, ch, af) as below\n")
pi <- sapply(trios, function(t) {
  d <- data.frame(L.1 = t[c(1, 3, 5)], L.2 = t[c(2, 4, 6)]); rownames(d) <- persons
  l <- lr(list(pf, pu), d); l[1] / l[2]
})
carr("ref_pi", pi)

# 2. kinship index: persons a, b; full sibs, half sibs, parent-child vs
# unrelated (theta 0)
fs <- FamiliasPedigree(id = c("a", "b", "f", "m"), dadid = c("f", "f", NA, NA), momid = c("m", "m", NA, NA),
                       sex = c("male", "male", "male", "female"))
hs <- FamiliasPedigree(id = c("a", "b", "f", "m1", "m2"), dadid = c("f", "f", NA, NA, NA),
                       momid = c("m1", "m2", NA, NA, NA), sex = c("male", "male", "male", "female", "female"))
po <- FamiliasPedigree(id = c("a", "b"), dadid = c(NA, "a"), momid = c(NA, NA), sex = c("male", "male"))
un <- FamiliasPedigree(id = c("a", "b"), dadid = c(NA, NA), momid = c(NA, NA), sex = c("male", "male"))
pairs <- list(c("A","B", "A","B"), c("A","A", "A","B"), c("A","B", "C","D"), c("C","D", "C","C"),
              c("E","E", "E","E"), c("B","F", "F","B"), c("A","C", "C","E"))
cat("// Familias: KI for full sibs, half sibs, parent-child against unrelated, pairs (a, b) as below\n")
ki <- t(sapply(pairs, function(g) {
  d <- data.frame(L.1 = g[c(1, 3)], L.2 = g[c(2, 4)]); rownames(d) <- c("a", "b")
  l <- lr(list(fs, hs, po, un), d); l[1:3] / l[4]
}))
carr("ref_ki", ki)   # column-major: all full-sib values, then half sibs, then parent-child

# 3. match probability with the coancestry correction: the joint
# likelihood of two unrelated people sharing genotype G under kinship
# theta, over the single-person likelihood under the same theta
one <- FamiliasPedigree(id = "a", dadid = NA, momid = NA, sex = "male")
gts <- list(c("A","A"), c("A","B"), c("C","F"), c("F","F"), c("D","E"))
th <- c(0, 0.01, 0.03, 0.1)
cat("// Familias: P(G, G | unrelated, kinship th) / P(G | th), genotypes x theta as below\n")
mp <- sapply(th, function(k) sapply(gts, function(g) {
  d2 <- data.frame(L.1 = g[c(1, 1)], L.2 = g[c(2, 2)]); rownames(d2) <- c("a", "b")
  d1 <- data.frame(L.1 = g[1], L.2 = g[2]); rownames(d1) <- "a"
  lr(list(un), d2, k) / lr(list(one), d1, k)
}))
carr("ref_mp", mp)   # genotypes fastest

# 4. the number of distinct alleles among 2x independent alleles at the
# locus, x = 1..3 contributors (theta 0)
cat("// numberofalleles: pr_number_of_distinct_alleles(2 x, fr), x = 1, 2, 3 (counts 1..6)\n")
for (x in 1:3) {
  v <- pr_number_of_distinct_alleles(2 * x, fr)
  carr(sprintf("ref_nd%d", x), as.numeric(v)[1:6])
}

# 5. Cllr and its minimum over monotone recalibrations: comparison's
# empirical cross-entropy at prior 0.5 (base-2 logs)
lss <- exp(seq(-1, 3, length.out = 15)); lds <- exp(seq(-4, 1.2, length.out = 18))
e <- calc.ece(lss, lds)
cat("// comparison: calc.ece(lss, lds) at prior 0.5: ece (Cllr) and ece.cal (Cllr of the PAV-calibrated LRs)\n")
carr("ref_ece", c(e$ece[e$prior == 0.5], e$ece.cal[e$prior == 0.5]))
cat("// lss = exp(seq(-1, 3, length.out = 15)), lds = exp(seq(-4, 1.2, length.out = 18))\n")

# 6. DNAtools: the same distribution by its compiled routine (theta 0)
cat("DNAtools", format(packageVersion("DNAtools")), "\n")
cat("// DNAtools: Pnm_locus(m = x, theta = 0, alleleProbs = fr), x = 1, 2, 3 (counts 1..2x)\n")
for (x in 1:3) carr(sprintf("ref_dna%d", x), Pnm_locus(m = x, theta = 0, alleleProbs = fr))

# 7. fairmodels (DALEX explainer around fixed scores): per-group rates of
# the deterministic sample of dev/doss/fairness.R, at cutoff 0.5 and at
# 0.001 (almost every case predicted positive; its cutoff must lie in [0, 1])
suppressMessages({library(DALEX); library(fairmodels)})
cat("fairmodels", format(packageVersion("fairmodels")), "\n")
i <- 0:999; y <- as.integer(((i * 13) %% 7) < 3); g <- c("a", "b", "c")[i %% 3 + 1]
p <- (((i * 37) %% 101) + 30 * y) / 131
e <- explain(list(), data = data.frame(x = p, y = y), y = y, predict_function = function(m, nd) nd$x,
             label = "s", verbose = FALSE, model_info = list(package = "none", ver = "0", type = "classification"))
for (ct in c(0.5, 0.001)) {
  fc <- fairness_check(e, protected = factor(g), privileged = "a", cutoff = ct, verbose = FALSE)
  gd <- fc$groups_data$s
  cat(sprintf("// fairness_check(e, protected, \"a\", cutoff = %g)$groups_data$s: ACC TPR TNR FPR FNR PPV NPV STP (groups a b c)\n", ct))
  carr(sprintf("ref_fm%s", if (ct == 0.5) "" else "_low"), unlist(lapply(c("ACC", "TPR", "TNR", "FPR", "FNR", "PPV", "NPV", "STP"), function(k) gd[[k]])))
}
