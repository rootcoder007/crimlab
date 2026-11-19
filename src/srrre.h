// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// rOpenSci statistical software standards for regression software, read
// by the srr package (srr::srr_report()), as src/srr.h holds the General
// and EDA ones and src/srrsp.h the spatial ones. They are about the fits
// of ppmlfe(), ivfe(), ppmlcf(), expgmm(), fepnl(), gmmpnl(), didfe(),
// rdloc() and sarml(). Standards met by a test are tagged there
// (tests/testthat/test-racc.R); the rest are here. Nothing here is
// compiled: no source file includes this header.
#ifndef CLAB_SRRRE_H
#define CLAB_SRRRE_H

//' srr_stats
//'
//' Regression standards the package meets.
//'
//' @srrstats {RE1.0} There is no formula interface, and the help page
//'   `racc` ("How a model is given") gives the reason: a model is given
//'   by the names of columns, each in the argument of its role (outcome,
//'   regressors, instruments, fixed effects, clusters), which a formula
//'   has no place for, and the fixed effects are absorbed, not expanded
//'   into columns.
//' @srrstats {RE1.2} Each help page says what its regressors are (the
//'   names of numeric columns of `data`, or a numeric vector or matrix),
//'   and `racc` says what is not taken: a column of labels is a fixed
//'   effect or a cluster, never a regressor.
//' @srrstats {RE1.3} The terms of a fit carry the names of the input
//'   columns into its table of estimates and its covariance matrix.
//' @srrstats {RE1.3a} `racc` ("What a fit keeps") says what is not
//'   carried: row names, since a fit has one row per term.
//' @srrstats {RE1.4} Each help page states the model it fits and what its
//'   standard errors allow (clustered, two-way, heteroskedastic,
//'   Newey-West, bootstrap); ppmlfe() documents separation, where the
//'   estimates do not exist, with its three ways of handling it, tested
//'   in test-apx.R; ivfe() reports the strength of the instruments it
//'   assumes relevant (first-stage F statistics with a critical value).
//' @srrstats {RE2.0} The one transformation is the absorbing of the fixed
//'   effects, documented on each help page; a model without them is
//'   another model, so nothing switches it off. Labels given as fixed
//'   effects or clusters are read as groups, as each help page says.
//' @srrstats {RE2.3} The count models take an offset, the exposure
//'   (`expo`): the mean is the exposure times the fitted rate. It is
//'   documented on their help pages and tested in test-apx.R and
//'   test-didfe.R.
//' @srrstats {RE2.4} Collinearity is looked for before any estimate is
//'   returned; see RE2.4a and RE2.4b in test-racc.R.
//' @srrstats {RE3.2} The iterative fits stop at a relative change of
//'   `tol` (1e-12 for ppmlfe() and ppmlcf(), 1e-10 for expgmm()) or
//'   after `maxit` iterations (100); both defaults are in the usage of
//'   each help page and explained there.
//' @srrstats {RE3.3} `tol` and `maxit` are arguments.
//' @srrstats {RE4.0} Each fit is a list with a class of its own, named
//'   after the function (`ppmlfe`, `ivfe`, `ppmlcf`, `expgmm`, `fepnl`,
//'   `gmmpnl`, `didfe`, `rdloc`, `sarml`).
//' @srrstats {RE4.7} The iterative fits return `conv` (whether the fit
//'   converged) and `iter` (the iterations used).
//' @srrstats {RE4.10} fepnl() and gmmpnl() return their residuals
//'   (`resid`), described on their help pages. The other fits do not:
//'   their residuals would need the absorbed fixed effects (RE4.9).
//' @srrstats {RE4.11} The fits return statistics of fit: the deviance and
//'   pseudo log-likelihood (ppmlfe()), the residual sum of squares and
//'   the first-stage and over-identification statistics (ivfe()), R
//'   squared and within R squared (fepnl()), Hansen's J (expgmm()), the
//'   Arellano-Bond and Sargan-Hansen tests (gmmpnl()), the log-likelihood
//'   and AIC (sarml()).
//' @srrstats {RE4.13} The regressors are named in the table of estimates
//'   (`coef$term`), and the fits record the fixed effects, the clusters
//'   and the kind of standard error used.
//' @srrstats {RE4.17} Each of the nine classes has a print method that
//'   shows what was fitted (rows, fixed effects, standard errors) and the
//'   table of estimates.
//' @srrstats {RE4.18} ppmlfe(), ivfe(), ppmlcf(), expgmm() and rdloc()
//'   have summary methods, which give the table of estimates; no summary
//'   needs a computation of its own.
//' @srrstats {RE6.0} Each of the nine classes has a plot method.
//' @srrstats {RE6.1} They are S3 methods, `plot.<class>`.
//' @noRd
#define CLAB_SRR_RE 1

//' NA_standards
//'
//' Regression standards that do not apply, with the reason for each.
//'
//' @srrstatsNA {RE1.1} No formula interface (RE1.0), so no conversion of
//'   a formula to a matrix to document.
//' @srrstatsNA {RE2.2} Missing values are refused in the outcome and in
//'   the regressors alike (RE2.1, G2.14b), and no fit predicts, so there
//'   is no case of fitting on complete predictors to fill in a missing
//'   response.
//' @srrstatsNA {RE4.1} No fit is so long that a model object without the
//'   fit would help; the standard is optional.
//' @srrstatsNA {RE4.4} No formula interface (RE1.0): the specification is
//'   in the fit's own fields (the terms, the fixed effects, the
//'   clusters), shown by print().
//' @srrstatsNA {RE4.8} The fits estimate coefficients of a model whose
//'   fixed effects are absorbed; they do not return the response.
//' @srrstatsNA {RE4.9} Modelled values need the fixed effects, which are
//'   absorbed and never estimated, so the fits do not return them; the
//'   help page `racc` says so.
//' @srrstatsNA {RE4.12} No fit transforms its input beyond absorbing the
//'   fixed effects, which has no inverse to offer.
//' @srrstatsNA {RE4.14} No fit forecasts or extrapolates.
//' @srrstatsNA {RE4.15} No fit forecasts.
//' @srrstatsNA {RE4.16} There is no predict method.
//' @srrstatsNA {RE6.2} The fits return no fitted values (RE4.9); the plot
//'   of a fit draws its estimates with their intervals.
//' @srrstatsNA {RE6.3} No fit forecasts.
//' @srrstatsNA {RE7.4} No fit forecasts.
//' @noRd
#endif
