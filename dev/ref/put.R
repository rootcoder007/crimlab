# Shared by the generators in dev/ref: writes the list x of reference
# values to tests/testthat/ref/<name>.txt as one R expression with hex
# floats (exact), headed by the packages and versions that computed it,
# and checks that it reads back identical. The tests read it with refget().
put <- function(x, name, pkgs) {
  v <- vapply(pkgs, function(p) as.character(utils::packageVersion(p)), "")
  x$pkgs <- v
  f <- file.path("tests", "testthat", "ref", paste0(name, ".txt"))
  dir.create(dirname(f), showWarnings = FALSE)
  co <- c("keepNA", "keepInteger", "niceNames", "showAttributes", "hexNumeric")
  writeLines(c(
    sprintf("# Written by dev/ref/%s.R with %s.", name,
            paste(pkgs, v, collapse = ", ")),
    "# Do not edit: run the generator again.",
    deparse(x, control = co)
  ), f)
  stopifnot(identical(dget(f), x))
  cat(f, "written:", length(unlist(x)), "values\n")
}
