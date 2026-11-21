#' @srrstats {G1.4} Every exported function is documented with roxygen2.
#' @section Notice of a newer release:
#' The package makes no request on its own. [updchk()] asks the release
#' page once, when called. A notice at attach is off unless turned on:
#' with `options(crimlab.update = TRUE)` before the package is attached
#' (in `.Rprofile`, say), or the environment variable
#' `CRIMLAB_UPDATE_CHECK` set to `1`, an interactive session then says in
#' one line when a newer release exists. The page of the latest release
#' is asked for at most once a day, with a limit of two seconds, and the
#' answer is kept in `tools::R_user_dir("crimlab", "cache")`; nothing is
#' sent but that request, and a failure is silent. A session that is not
#' interactive (a script, a check, a test) never asks.
#' `options(crimlab.update = FALSE)`, or unsetting the variable, turns the
#' notice off again.
#' @keywords internal
#' @useDynLib crimlab, .registration = TRUE
"_PACKAGE"
