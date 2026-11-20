#' @srrstats {G1.4} Every exported function is documented with roxygen2.
#' @section Notice of a newer release:
#' In an interactive session, attaching the package says in one line when
#' a newer release exists. The page of the latest release is asked for at
#' most once a day, with a limit of two seconds, and the answer is kept in
#' `tools::R_user_dir("crimlab", "cache")`; nothing is sent but that
#' request, and a failure is silent. A session that is not interactive
#' (a script, a check, a test) never asks. `options(crimlab.update = FALSE)`
#' before the package is attached, or the environment variable
#' `CRIMLAB_NO_UPDATE_CHECK` set to anything, turns it off.
#' @keywords internal
#' @useDynLib crimlab, .registration = TRUE
"_PACKAGE"
