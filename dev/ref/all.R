# Runs every generator in dev/ref, each in its own R process, writing the
# fixtures to OUTDIR instead of tests/testthat/ref. A generator whose
# packages (its `need`) are not installed is skipped, with a line saying
# so. From the package root:
#   Rscript dev/ref/all.R OUTDIR
# Writes OUTDIR/<name>.txt, OUTDIR/<name>.log and OUTDIR/status.csv (name,
# status ok, skipped or failed, detail); exits 1 if a generator failed.
a <- commandArgs(TRUE)
if (length(a) != 1L) stop("usage: Rscript dev/ref/all.R OUTDIR")
dir.create(a[1], showWarnings = FALSE, recursive = TRUE)
out <- normalizePath(a[1])
gen <- sub("\\.R$", "", list.files(file.path("dev", "ref"), "\\.R$"))
gen <- setdiff(gen, c("all", "cmp", "put"))
rs <- file.path(R.home("bin"), "Rscript")
Sys.setenv(REFOUT = out)
st <- data.frame(name = gen, status = "", detail = "")
for (i in seq_along(gen)) {
  f <- file.path("dev", "ref", paste0(gen[i], ".R"))
  need <- eval(str2lang(sub("^need <- ", "", grep("^need <- ", readLines(f),
                                                  value = TRUE)[1])))
  ok <- vapply(need, requireNamespace, NA, quietly = TRUE)
  if (!all(ok)) {
    st$status[i] <- "skipped"
    st$detail[i] <- paste("not installed:", paste(need[!ok], collapse = " "))
  } else {
    log <- file.path(out, paste0(gen[i], ".log"))
    rc <- system2(rs, f, stdout = log, stderr = log)
    st$status[i] <- if (identical(rc, 0L)) "ok" else "failed"
    st$detail[i] <- paste(need, vapply(need, function(p) {
      as.character(utils::packageVersion(p))
    }, ""), collapse = ", ")
  }
  message(sprintf("%-8s %-8s %s", gen[i], st$status[i], st$detail[i]))
}
utils::write.csv(st, file.path(out, "status.csv"), row.names = FALSE)
quit(status = as.integer(any(st$status == "failed")))
