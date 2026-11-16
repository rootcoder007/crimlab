#' Rao-Scott chi-square test for a survey table
#'
#' The test of independence of two categorical variables in a two-way
#' table from a stratified, clustered and weighted sample (a victimisation
#' or justice survey, say), with Rao and Scott's corrections for the
#' design, to go with the estimates of [svymx()] and the tables of
#' [ctabx()]. Pearson's statistic on such data is too large: the design
#' makes the cell proportions vary more than simple random sampling, and
#' the corrections divide it by the design effects.
#'
#' The cell proportions `p` are weighted; Pearson's `X^2` is computed on
#' them scaled to the number of respondents `n`, `X^2 = n sum (p_ij -
#' p_i. p_.j)^2 / (p_i. p_.j)`. Their covariance `V` comes by linearisation
#' over the first-stage units within strata, as [svymx()] reads the design
#' (Lumley 2004). The generalised design effects `deff` are the
#' eigenvalues of `Delta = (C' D^-1 C / n)^-1 (C' D^-1 V D^-1 C)`, `C` a
#' basis of the interaction contrasts and `D` the diagonal of `p` (Rao and
#' Scott 1984). The tests: `pearson`, `X^2` on `K = (r - 1)(c - 1)`
#' degrees of freedom, wrong under the design; `rs1`, the first-order
#' correction `X^2 / mean(deff)` on `K`; `rs2`, the second-order
#' (Satterthwaite) correction `X^2 / (mean(deff) (1 + a^2))` on `K / (1 +
#' a^2)`, `a` the coefficient of variation of the design effects; and `f`,
#' `X^2 / sum(deff)` as an F on `d0 = K / (1 + a^2)` and `d0 nu`
#' degrees of freedom, `nu` the first-stage units less the strata, the
#' form Thomas and Rao (1987) found to hold its level best. They agree
#' with `survey::svychisq()` (`"Chisq"` for `rs1`, `"F"` for `f`).
#'
#' @param x,y The two variables: labels (numbers, text or a factor), one
#'   per respondent, each with two labels or more; for `print()`, `x` is an
#'   object of class `rscchi`.
#' @inheritParams svymx
#' @return An object of class `rscchi`, a list: `tab`, the weighted cell
#'   proportions (a matrix, `x` by `y`); `test`, a data frame (`test`:
#'   `pearson`, `rs1`, `rs2`, `f`; `stat`, `df1`, `df2`, `p`); `deff`, the
#'   generalised design effects; `dbar` and `a2`, their mean and squared
#'   coefficient of variation; and `des` (`n`, `strata`, `psu`, `df`).
#' @references
#'   Rao, J. N. K. and Scott, A. J. (1981). The analysis of categorical data
#'   from complex sample surveys: chi-squared tests for goodness of fit and
#'   independence in two-way tables. *Journal of the American Statistical
#'   Association*, 76(374), 221-230. \doi{10.1080/01621459.1981.10477633}
#'
#'   Rao, J. N. K. and Scott, A. J. (1984). On chi-squared tests for
#'   multiway contingency tables with cell proportions estimated from survey
#'   data. *The Annals of Statistics*, 12(1), 46-60.
#'   \doi{10.1214/aos/1176346391}
#'
#'   Thomas, D. R. and Rao, J. N. K. (1987). Small-sample comparisons of
#'   level and power for simple goodness-of-fit statistics under cluster
#'   sampling. *Journal of the American Statistical Association*, 82(398),
#'   630-636. \doi{10.1080/01621459.1987.10478476}
#'
#'   Lumley, T. (2004). Analysis of complex survey samples. *Journal of
#'   Statistical Software*, 9(8), 1-19. \doi{10.18637/jss.v009.i08}
#' @examples
#' set.seed(1)
#' d <- data.frame(region = rep(c("north", "south"), each = 60),
#'                 block = rep(1:12, each = 10),
#'                 age = sample(c("young", "middle", "old"), 120, TRUE),
#'                 w = runif(120, 50, 150))
#' d$victim <- rbinom(120, 1, ifelse(d$age == "young", 0.35, 0.2))
#' rscchi(d$age, d$victim, d$w, d$region, d$block)
#' @export
rscchi <- function(x, y, w = NULL, strata = NULL, psu = NULL, fpc = NULL) {
  k <- .Call(C_clab_r_rscchi, x, y, w, strata, psu, fpc)
  nr <- length(k$xlev)
  tab <- matrix(k$tab, nr, dimnames = list(k$xlev, k$ylev))
  t <- k$test
  t[is.nan(t)] <- NA_real_
  a <- k$aux
  structure(
    list(tab = tab,
         test = data.frame(test = c("pearson", "rs1", "rs2", "f"),
                           stat = t[1:4], df1 = t[5:8], df2 = t[9:12],
                           p = t[13:16], stringsAsFactors = FALSE),
         deff = k$deff, dbar = a[1L], a2 = a[2L],
         des = c(n = length(x), strata = k$nstrata, psu = k$npsu,
                 df = a[3L])),
    class = "rscchi"
  )
}

#' @rdname rscchi
#' @param ... Not used by `print()`.
#' @export
print.rscchi <- function(x, ...) {
  g <- x[["des"]]
  cat("rscchi: ", nrow(x[["tab"]]), " x ", ncol(x[["tab"]]), " table, ",
      g[["n"]], " respondents, ", g[["strata"]], " strata, ", g[["psu"]],
      " first-stage units\n", sep = "")
  print(x[["test"]], row.names = FALSE, digits = 4)
  invisible(x)
}
