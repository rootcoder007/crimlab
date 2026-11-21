# Name gate: stops when an exported name is already a function in a CRAN
# package. Run from the package root before adding an export.
ex <- getNamespaceExports(asNamespace("crimlab"))
hit <- get("fn_names_on_cran", asNamespace("pkgcheck"))(ex)
# rmoriebricklayer is by the same author, and the functions crimlab took
# over from it keep their names
hit <- hit[hit[["package"]] != "rmoriebricklayer", , drop = FALSE]
if (nrow(hit)) {
  print(hit)
  stop("exported names taken on CRAN", call. = FALSE)
}
cat(length(ex), "exported names free on CRAN\n")
