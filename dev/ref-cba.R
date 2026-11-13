# The stored reference values of tests/testthat/test-cbatt.R, from
# 'WeightIt', 'MatchIt', 'cobalt' and 'EValue' on the data of
# tests/testthat/helper-cba.R. From the package root:
#   Rscript dev/ref-cba.R
# The tests hold the printed numbers; these packages are not needed to
# run them.
sys.source(file.path("tests", "testthat", "helper-cba.R"), envir = environment())
z <- cbdata()
dd <- cbind(z$x, d = z$d, y = z$y)
p <- function(nm, x) cat(nm, "=", deparse(signif(unname(x), 15)), "\n")
for (pk in c("WeightIt", "MatchIt", "cobalt", "EValue")) {
  cat(pk, format(packageVersion(pk)), "\n")
}
dm <- function(w) {
  weighted.mean(z$y[z$d == 1], w[z$d == 1]) -
    weighted.mean(z$y[z$d == 0], w[z$d == 0])
}
for (e in c("ATE", "ATT", "ATC")) {
  w <- WeightIt::weightit(d ~ a + b, data = dd, method = "glm", estimand = e)
  p(paste0("ipw.", e), dm(w$weights))
}
eb <- WeightIt::weightit(d ~ a + b, data = dd, method = "ebal", estimand = "ATT")
p("ebal", dm(eb$weights))
m <- MatchIt::matchit(d ~ a + b, data = dd, method = "nearest", caliper = 0.2)
md <- MatchIt::match.data(m)
p("psm", weighted.mean(md$y[md$d == 1], md$weights[md$d == 1]) -
  weighted.mean(md$y[md$d == 0], md$weights[md$d == 0]))
p("psm.n", nrow(md))
w <- WeightIt::weightit(d ~ a + b, data = dd, method = "glm", estimand = "ATE")
b <- cobalt::bal.tab(d ~ a + b, data = dd, weights = w$weights,
                     s.d.denom = "pooled", disp.v.ratio = TRUE, un = TRUE)
p("smd", b$Balance$Diff.Un)
p("vr", b$Balance$V.Ratio.Un)
p("smdw", b$Balance$Diff.Adj)
p("vrw", b$Balance$V.Ratio.Adj)
p("ev.rr", EValue::evalues.RR(3.9, 1.8, 8.7)[2, 1:2])
p("ev.hr", EValue::evalues.HR(0.8, 0.71, 0.91, rare = FALSE)[2, c(1, 3)])
p("ev.or", EValue::evalues.OR(2.5, 1.5, 4, rare = FALSE)[2, 1:2])
p("ev.md", EValue::evalues.MD(0.5, 0.1)[2, 1:2])
