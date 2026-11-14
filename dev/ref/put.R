# Shared by the generators in dev/ref: writes the list x of reference
# values to <dir>/<name>.txt as one R expression with hex floats (exact),
# headed by the packages and versions that computed it and the platform
# they ran on, and checks that it reads back identical. dir is REFOUT
# when set (dev/ref/all.R sets it), tests/testthat/ref otherwise. The
# tests read the file with refget().
put <- function(x, name, pkgs) {
  v <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
  x$pkgs <- v
  out <- Sys.getenv("REFOUT", file.path("tests", "testthat", "ref"))
  f <- file.path(out, paste0(name, ".txt"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  co <- c("keepNA", "keepInteger", "niceNames", "showAttributes", "hexNumeric")
  writeLines(c(
    sprintf("# Written by dev/ref/%s.R with %s.", name,
            paste(pkgs, v, collapse = ", ")),
    sprintf("# On %s, %s, %s.", R.version.string, R.version$platform,
            utils::sessionInfo()$running),
    "# Do not edit: run the generator again.",
    deparse(x, control = co)
  ), f)
  stopifnot(identical(dget(f), x))
  cat(f, "written:", length(unlist(x)), "values\n")
}
