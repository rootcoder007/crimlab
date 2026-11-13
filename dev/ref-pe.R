# The stored reference values of the policy-evaluation tests
# (tests/testthat/test-didfe.R, test-evtst.R, test-stgdid.R), from
# 'fixest', 'did' and 'bacondecomp' on the panel of
# tests/testthat/helper-pe.R. From the package root:
#   Rscript dev/ref-pe.R
# The tests hold the printed numbers; these packages are not needed to
# run them.
sys.source(file.path("tests", "testthat", "helper-pe.R"), envir = environment())
d <- pepanel()
p <- function(nm, x) cat(nm, "=", deparse(signif(unname(x), 15)), "\n")
cat("fixest", format(packageVersion("fixest")), "\n")
f1 <- fixest::feols(rate ~ on + x | town + year, d, cluster = ~region)
p("ols1", fixest::se(f1))
f2 <- fixest::feols(rate ~ on + x | town + year, d, cluster = ~ region + year)
p("ols2", fixest::se(f2))
f3 <- fixest::feols(rate ~ on | town + year, d, weights = ~w, cluster = ~town)
p("olsw", c(coef(f3), fixest::se(f3)))
f4 <- fixest::fepois(cnt ~ on | town + year, d, offset = ~log(pop),
                     cluster = ~town)
p("pois", c(coef(f4), fixest::se(f4)))
f5 <- fixest::fepois(cnt ~ on | town + year, d, cluster = ~ region + year)
p("pois2", c(coef(f5), fixest::se(f5)))
d$k <- ifelse(is.na(d$first), -1000, pmin(pmax(d$year - d$first, -3), 3))
f6 <- fixest::feols(rate ~ i(k, ref = c(-1, -1000)) | town + year, d,
                    cluster = ~town)
p("evt", c(coef(f6), fixest::se(f6)))
w6 <- fixest::wald(f6, keep = "k::-3|k::-2", print = FALSE)
p("evtf", c(w6$stat, w6$p))
cat("bacondecomp", format(packageVersion("bacondecomp")), "\n")
b <- bacondecomp::bacon(rate ~ on, d, id_var = "town", time_var = "year",
                        quietly = TRUE)
b <- b[order(b$type, b$treated, b$untreated), ]
print(b)
cat("did", format(packageVersion("did")), "\n")
d$g <- ifelse(is.na(d$first), 0, d$first)
for (ctl in c("notyettreated", "nevertreated")) {
  a <- did::att_gt("rate", "year", "town", "g", data = d,
                   control_group = ctl, bstrap = FALSE, cband = FALSE)
  p(paste0("att.", ctl), a$att)
}
for (m in c("dr", "ipw", "reg")) {
  a <- did::att_gt("rate", "year", "town", "g", xformla = ~x, data = d,
                   control_group = "notyettreated", est_method = m,
                   bstrap = FALSE, cband = FALSE)
  p(paste0("attx.", m), a$att)
}
