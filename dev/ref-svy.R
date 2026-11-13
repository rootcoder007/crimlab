# The stored reference values of tests/testthat/test-svymx.R, from the
# 'survey' package on the data of tests/testthat/helper-svy.R. From the
# package root:
#   Rscript dev/ref-svy.R
sys.source(file.path("tests", "testthat", "helper-svy.R"), envir = environment())
d <- svdata()
p <- function(nm, x) cat(nm, "=", deparse(signif(unname(x), 15)), "\n")
cat("survey", format(packageVersion("survey")), "\n")
ds <- survey::svydesign(ids = ~block, strata = ~region, weights = ~w, data = d)
m <- survey::svymean(~y, ds)
p("mean", c(coef(m), survey::SE(m)))
tt <- survey::svytotal(~y, ds)
p("total", c(coef(tt), survey::SE(tt)))
r <- survey::svyratio(~y, ~x, ds)
p("ratio", c(coef(r), survey::SE(r)))
v <- survey::svyvar(~y, ds)
p("var", c(v, survey::SE(v)))
pr <- survey::svymean(~v, ds)
p("prop", c(coef(pr), survey::SE(pr)))
b <- survey::svyby(~y, ~region, ds, survey::svymean, covmat = TRUE)
p("by", c(coef(b), survey::SE(b)))
k <- survey::svycontrast(b, c(-1, 1, 0))
p("cmp", c(coef(k), survey::SE(k)))
dsf <- survey::svydesign(ids = ~block, strata = ~region, weights = ~w,
                         fpc = ~N, data = d)
p("fpc", survey::SE(survey::svymean(~y, dsf)))
