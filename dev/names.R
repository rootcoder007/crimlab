# Name gate: stops when an exported name is already a function in a CRAN
# package. Run from the package root before adding an export.
ex <- getNamespaceExports(asNamespace("crimlab"))
hit <- get("fn_names_on_cran", asNamespace("pkgcheck"))(ex)
if (nrow(hit)) {
  print(hit)
  stop("exported names taken on CRAN", call. = FALSE)
}
cat(length(ex), "exported names free on CRAN\n")
