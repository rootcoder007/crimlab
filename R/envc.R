#' LMDI decomposition of a count over the factors of a Kaya identity
#'
#' How much of a change in recorded environmental offences (pollution
#' cases, waste violations) between two years, or of a region's gap to a
#' reference region, is owed to each factor of an identity that multiplies
#' out to the count: the additive log-mean Divisia index (LMDI)
#' decomposition, temporal and spatial, as Duan et al. (2026, eqs 1-5,
#' pp. 8-9) apply it to crimes of polluting the environment in 30 Chinese
#' provinces, 2016-2023.
#'
#' The identity is built as a chain from the quantities `q_1, ..., q_m`
#' (the columns of `q`, in order): `y = (y / q_1) (q_1 / q_2) ... (q_{m-1}
#' / q_m) q_m`, `m + 1` factors that multiply to `y` whatever the
#' quantities. Duan et al.'s eq 1 is the chain of population `P`, area
#' `S`, emissions `W`, GDP, population again and the unemployed `UN`:
#' `CPEs = CPER x PD x EIPA^-1 x PEI x EDL x UR^-1 x UN`, the case rate
#' per head, population density, area per unit of emission, emission per
#' unit of GDP, GDP per head, people per unemployed person and the
#' unemployed (see the example). The factor values themselves, the rates
#' of the identity, are returned in `fac`.
#'
#' With `type = "time"` the change of each unit between consecutive
#' periods of `at` is `y_T - y_O = sum_j L(y_T, y_O) ln(f_jT / f_jO)` (eqs
#' 2-3), `L(a, b) = (a - b) / (ln a - ln b)` the logarithmic mean; a unit
#' missing either period is left out. With `type = "space"` each unit of
#' each period in `at` is set against a reference region of the same
#' period (eqs 4-5): with `ref = "mean"`, the region whose count and
#' quantities are the means over the units present, so its count is the
#' average count Duan et al. benchmark against; with `ref = "gmean"`, the
#' geometric means. Duan et al. do not say how the reference region's
#' factors were formed; the arithmetic means of the quantities are this
#' function's reading of "the national average", and `"gmean"` the
#' alternative. The decomposition has no residual: the contributions add
#' up to `dy`.
#'
#' A zero count or quantity has no logarithm. With `zero = 0` (Duan et
#' al.'s "Null: no numerical value", figs 5-6) the rows it touches get
#' `NA` contributions; with `zero` positive, zeros are replaced by it
#' first. The contributions then depend on `zero`, and only slowly fade as
#' it shrinks (as one over its logarithm), so a replaced zero should be
#' reported as such.
#'
#' @section Environmental crime: The family treats offences against the
#'   environment as the outcome: where, when and why pollution offences,
#'   illegal waste disposal and trafficking, and safety crimes in
#'   extractive industries are recorded, and how enforcement behaves. Its
#'   tools are [lmdix()] (decomposition over a Kaya identity, Duan et al.
#'   2026), [fepnl()] (two-way fixed effects with an interaction and the
#'   checks Duan et al. report) and [pcdtst()] (cross-sectional dependence
#'   of a panel, Andreatta et al. 2023), with the package's spatial tests
#'   ([gmoran()], [gearc()], [lmorn()]) for the spatial checks Andreatta et
#'   al. run. Recorded counts are counts of detected, prosecuted or judged
#'   cases: Duan et al. count first-instance judgments, Andreatta et al.
#'   the violations verified by law enforcement and gathered by an NGO, and
#'   both say the counts carry enforcement as well as offending (Andreatta
#'   et al. write the observed rate as the true rate plus a measurement
#'   error that can move with enforcement, p. 1372). Safety crimes are
#'   "violations of law by employers that either do, or have the potential
#'   to, cause sudden death or injury as a result of work-related
#'   activities" (Tombs and Whyte, quoted by Pitoulis 2014, p. 2); they are
#'   mostly regulatory offences enforced by inspectorates rather than
#'   police, mostly end in guilty pleas and fines (p. 47), and the injuries
#'   behind them are under-reported (p. 3, note 1), so recorded counts
#'   understate them. None of the three sources specifies a model of
#'   recorded against true offences, and the family has none.
#' @param y The count of each row (a unit in a period): numbers from 0.
#' @param q The quantities of the chain, in order: a numeric matrix or
#'   data frame with one row per row of `y` and one column per quantity,
#'   numbers from 0; the column names name the factors.
#' @param unit The unit of each row: labels.
#' @param time The period of each row: numbers; each unit at most once
#'   per period.
#' @param type `"time"` (change between periods) or `"space"` (gap to the
#'   reference region of each period).
#' @param at `NULL` or the periods to use, increasing: for `"time"`, the
#'   end points of consecutive intervals (default the first and last
#'   period; Duan et al. used 2016, 2018, 2020 and 2023); for `"space"`,
#'   the periods to decompose (default all).
#' @param ref `"mean"` or `"gmean"`: the reference region of `"space"`.
#' @param zero `0` (rows with a zero get `NA`) or a positive number put in
#'   place of the zeros.
#' @param fnm `NULL` or the names of the `ncol(q) + 1` factors.
#' @return An object of class `lmdix`, a list: `dec`, one row per unit and
#'   interval (`unit`, `from`, `to`, `y0`, `y1`, `dy`) or per unit and
#'   period (`unit`, `time`, `ref` the reference count, `y`, `dy`),
#'   followed by one column of contributions per factor; `fac`, the count
#'   and the factor values of every row (`unit`, `time`, `y`, one column
#'   per factor); `ref`, for `"space"`, the reference region's count and
#'   factors by period, else `NULL`; `factors`, the factor names; `type`;
#'   and `refkind`.
#' @references
#'   Duan, H., Nie, C., He, B., Chen, X., Duan, Z., Chen, S., Song, J. and
#'   Gao, Z. (2026). A multi-method framework for decoding the
#'   spatiotemporal patterns and driving factors of environmental
#'   pollution crimes in China. *Humanities and Social Sciences
#'   Communications*. \doi{10.1057/s41599-026-08766-9}
#'
#'   Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
#'   Digging into waste: an analysis of waste crime in the Italian
#'   provinces. *Regional Studies*, 57(7), 1367-1379.
#'   \doi{10.1080/00343404.2022.2135697}
#'
#'   Pitoulis, T. T. (2014). *Bloody Oil: A Critical Discourse Analysis of
#'   Safety Crimes in the Alberta Oil and Gas Industry*. MA thesis in
#'   criminology, University of Ottawa.
#' @examples
#' # six provinces over four years: cases, population (millions), area
#' # (1000 km2), emissions (1000 t), GDP (billions) and the unemployed
#' # (millions), simulated
#' set.seed(1)
#' d <- expand.grid(year = c(2016, 2018, 2020, 2023), prov = LETTERS[1:6])
#' i <- match(d$prov, LETTERS)
#' d$P <- 10 * i * (1 + 0.01 * (d$year - 2016))
#' d$S <- 50 + 20 * i
#' d$W <- 30 * i * (1 - 0.03 * (d$year - 2016))
#' d$GDP <- 40 * i * (1 + 0.05 * (d$year - 2016))
#' d$UN <- 0.3 * i * (1 + 0.02 * (d$year - 2016))
#' d$cases <- rpois(24, 2 * i + 8 * (d$year == 2018))
#' q <- d[c("P", "S", "W", "GDP", "P", "UN")]
#' f <- c("CPER", "PD", "EIPA^-1", "PEI", "EDL", "UR^-1", "UN")
#' r <- lmdix(d$cases, q, d$prov, d$year, at = c(2016, 2018, 2020, 2023),
#'            fnm = f)
#' r
#' # the gap of each province to the average province, by year
#' lmdix(d$cases, q, d$prov, d$year, type = "space", fnm = f)$dec
#' @export
lmdix <- function(y, q, unit, time, type = "time", at = NULL, ref = "mean",
                  zero = 0, fnm = NULL) {
  if (is.data.frame(q)) q <- as.matrix(q)
  .Call(C_clab_r_lmdix, y, q, unit, time, type, at, ref, zero, fnm)
}

#' @export
print.lmdix <- function(x, ...) {
  d <- x[["dec"]]
  f <- x[["factors"]]
  cat("lmdix: ", x[["type"]], " decomposition, ", nrow(d), " rows, ",
      length(f), " factors\n", sep = "")
  tot <- colSums(d[c("dy", f)], na.rm = TRUE)
  cat("contributions summed over the rows decomposed:\n")
  print(tot)
  invisible(x)
}

#' @rdname lmdix
#' @param x An `lmdix` result.
#' @param ... Passed to the base-graphics call.
#' @param main,col,pal The title, colour and palette, as in [plotm].
#' @param units `NULL` (every unit) or the units whose rows are summed.
#' @details The plot draws the contribution of each factor summed over the
#'   rows (of `units`), as horizontal bars, and returns those sums.
#' @export
plot.lmdix <- function(x, ..., main = NULL, col = NULL, pal = NULL,
                       units = NULL) {
  .pck(x, "lmdix", c("dec", "factors"))
  d <- x$dec
  if (!is.null(units)) d <- d[d$unit %in% units, , drop = FALSE]
  v <- colSums(d[x$factors], na.rm = TRUE)
  .pbars(unname(v), names(v), .pmain(main, "LMDI contributions"),
         "Contribution to the change in the count", col, pal, ref = 0, ...)
}

#' Two-way fixed effects with an interaction, for counts of offences
#'
#' The regression Duan et al. (2026, eqs 6-8, pp. 9-11) use to check the
#' drivers of recorded pollution crimes: the count on continuous drivers
#' and the product of two of them, with unit and period fixed effects,
#' errors clustered by unit, standardised coefficients, variance inflation
#' factors, the simple slopes of one driver at values of the other (p.
#' 23), and the Jarque-Bera test of the residuals (p. 25).
#'
#' The model is `y_it = sum_j b_j x_itj + a_i + m_t + e_it`, fitted by least
#' squares after the within transformation (the kernel of [didfe()], so an
#' unbalanced panel is fine); with `int = c("a", "b")` the column `a:b`,
#' the product of `a` and `b`, is added to `x`. The variance is CR1 by
#' `clus` (`G / (G - 1) (n - 1) / (n - K)` times the sandwich, as in
#' [didfe()]); `ssc = "none"` drops the factor. Duan et al. fitted the
#' model with the 'linearmodels' PanelOLS estimator clustered by province
#' and do not state its small-sample factor. Tests and intervals use
#' Student's t on `G - 1` degrees of freedom.
#'
#' `beta` is the standardised coefficient `b SD_x / SD_y` of Duan et al.'s
#' step 4 (p. 11), with the standard deviations of the columns as given
#' (the interaction column for its term). `vif` is the variance
#' inflation factor `1 / (1 - R2_j)` of each column of `x` (with the
#' interaction) regressed on a constant and the others over the rows used;
#' Duan et al. drop drivers with a VIF above 10 before the fit (p. 10) and
#' do not say whether their VIFs were pooled, as here, or within. The
#' simple slope of `int[1]` at a value `m` of `int[2]` is `b_a + b_ab m`
#' with variance `V_aa + m^2 V_ab,ab + 2 m V_a,ab`; by default at the mean
#' and the mean plus and minus one standard deviation of `int[2]`, Duan et
#' al.'s points. `r2` is `1 - RSS / TSS` of the full model (unit and period
#' effects counted as fitted), `r2w` the within R-squared. The Jarque-Bera
#' statistic is `n / 6 (S^2 + (K - 3)^2 / 4)` of the residuals, with
#' moments divided by `n`, on two degrees of freedom.
#'
#' Duan et al.'s model is linear in a count. For a Poisson or negative
#' binomial model with the same effects see [didfe()] (binary treatment)
#' and [cntreg()]; Andreatta et al. (2023) fit a dynamic panel of the
#' logged rate by system GMM, which the package does not have.
#'
#' @inheritParams lmdix
#' @param y The outcome of each row: finite numbers.
#' @param x The drivers: a numeric matrix or data frame, one row per row
#'   of `y` (or one vector); column names name the terms. For `plot()`, a
#'   `fepnl` result.
#' @param int `NULL` or the names of two columns of `x` whose product is
#'   added: the slopes are those of the first at values of the second.
#' @param clus `NULL` (cluster by unit) or the cluster of each row
#'   (labels).
#' @param at `NULL` or the values of `int[2]` for the simple slopes.
#' @param ssc `"cr1"` or `"none"`: the small-sample factor of the
#'   variance.
#' @param level The confidence level of the intervals.
#' @return An object of class `fepnl`, a list: `coef`, one row per term
#'   (`term`, `est`, `se`, `lo`, `hi`, `t`, `p`, `beta`, `vif`); `slopes`
#'   (`at`, `slope`, `se`, `lo`, `hi`, `t`, `p`), `NULL` without `int`;
#'   `fit`, one row: `n`, `units`, `periods`, `ncl` (`G`), `k` (`K`), `r2`,
#'   `r2w`, `rss`; `jb`, a named vector `skew`, `kurt`, `jb`, `p`; `vcov`;
#'   `resid`, the residuals in the order of the rows; `ssc`; and `level`.
#' @references
#'   Duan, H., Nie, C., He, B., Chen, X., Duan, Z., Chen, S., Song, J. and
#'   Gao, Z. (2026). A multi-method framework for decoding the
#'   spatiotemporal patterns and driving factors of environmental
#'   pollution crimes in China. *Humanities and Social Sciences
#'   Communications*. \doi{10.1057/s41599-026-08766-9}
#'
#'   Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
#'   Digging into waste: an analysis of waste crime in the Italian
#'   provinces. *Regional Studies*, 57(7), 1367-1379.
#'   \doi{10.1080/00343404.2022.2135697}
#' @examples
#' # 30 provinces over 8 years: cases on the case rate, the unemployment
#' # rate and their product, simulated
#' set.seed(2)
#' d <- expand.grid(year = 2016:2023, prov = paste0("p", 1:30))
#' d$cper <- runif(240, 0.005, 0.04)
#' d$ur <- runif(240, 0.01, 0.03)
#' d$eipa <- rnorm(240, 5)
#' d$cases <- 50 + 9000 * d$cper + 1500 * d$ur - 4e5 * d$cper * d$ur +
#'   rnorm(240, sd = 5)
#' r <- fepnl(d$cases, d[c("cper", "eipa", "ur")], d$prov, d$year,
#'            int = c("cper", "ur"))
#' r
#' r$slopes
#' @export
fepnl <- function(y, x, unit, time, int = NULL, clus = NULL, at = NULL,
                  ssc = "cr1", level = 0.95) {
  if (is.data.frame(x)) x <- as.matrix(x)
  .Call(C_clab_r_fepnl, y, x, unit, time, int, clus, at, ssc, level)
}

#' @export
print.fepnl <- function(x, ...) {
  f <- x[["fit"]]
  cat("fepnl: ", f[["n"]], " rows, ", f[["units"]], " units, ",
      f[["periods"]], " periods, ", f[["ncl"]], " clusters, R2 ",
      format(f[["r2"]], digits = 4), ", within R2 ",
      format(f[["r2w"]], digits = 4), "\n", sep = "")
  print(x[["coef"]], row.names = FALSE)
  if (!is.null(x[["slopes"]])) {
    cat("simple slopes:\n")
    print(x[["slopes"]], row.names = FALSE)
  }
  j <- x[["jb"]]
  cat("Jarque-Bera of the residuals ", format(j[["jb"]], digits = 4),
      ", p ", format(j[["p"]], digits = 4), "\n", sep = "")
  invisible(x)
}

#' @rdname fepnl
#' @param ... Passed to the base-graphics call.
#' @param main,col,pal The title, colour and palette, as in [plotm].
#' @details The plot draws the coefficients with their intervals and, with
#'   an interaction, the simple slopes at each value of the moderator; it
#'   returns what it drew.
#' @export
plot.fepnl <- function(x, ..., main = NULL, col = NULL, pal = NULL) {
  .pck(x, "fepnl", "coef")
  op <- .ppanels(if (is.null(x$slopes)) 1L else 2L)
  on.exit(graphics::par(op), add = TRUE)
  r <- x$coef
  out <- list(coef = .pdots(r$term, r$est, r$lo, r$hi,
                            .pmain(main, "Two-way fixed-effects estimates"),
                            "Coefficient", col, pal, ref = 0, ...))
  s <- x$slopes
  if (!is.null(s)) {
    out$slopes <- .pdots(paste("at", format(s$at, digits = 3)), s$slope,
                         s$lo, s$hi, "Simple slopes", "Slope", col, pal,
                         ref = 0, ...)
  }
  invisible(out)
}

#' Pesaran's CD test of cross-sectional dependence in a panel
#'
#' Whether the residuals (or any variable) of a unit-by-period panel move
#' together across units, the test Andreatta et al. (2023, p. 1376) run on
#' the residuals of their waste-crime models before deciding a spatial
#' model is not needed. For each pair of units the correlation `rho_ab` of
#' their values over the `T_ab` periods they share, about the means over
#' those periods; `CD = sum sqrt(T_ab) rho_ab / sqrt(m)` over the `m`
#' pairs, which for a balanced panel of `N` units and `T` periods is
#' `sqrt(2 T / (N (N - 1))) sum rho_ab`, standard normal under no
#' dependence. Pairs that share fewer than two periods, or along which a
#' unit does not vary, are left out. This is the statistic of
#' `plm::pcdtest(test = "cd")`. For dependence between neighbours rather
#' than across all pairs, use [gmoran()] on each period's residuals, as
#' Andreatta et al. also do.
#'
#' @param x The values of each row: finite numbers (residuals of a panel
#'   model, typically).
#' @inheritParams lmdix
#' @return A named numeric vector: `cd`, `p` (two-sided normal), `rho` and
#'   `absrho` (the mean correlation and the mean absolute correlation of
#'   the pairs), `pairs` (`m`), `units` and `periods`.
#' @references
#'   Andreatta, D., Favarin, S., Lisciandra, M. and Millemaci, E. (2023).
#'   Digging into waste: an analysis of waste crime in the Italian
#'   provinces. *Regional Studies*, 57(7), 1367-1379.
#'   \doi{10.1080/00343404.2022.2135697}
#' @examples
#' set.seed(3)
#' d <- expand.grid(year = 2009:2017, prov = paste0("p", 1:20))
#' common <- rnorm(9)
#' d$e <- rnorm(180) + 0.8 * common[d$year - 2008]
#' pcdtst(d$e, d$prov, d$year)
#' # the residuals of a two-way fit
#' d$x <- rnorm(180)
#' d$y <- d$x + d$e
#' pcdtst(fepnl(d$y, d$x, d$prov, d$year)$resid, d$prov, d$year)
#' @export
pcdtst <- function(x, unit, time) {
  .Call(C_clab_r_pcdtst, x, unit, time)
}
