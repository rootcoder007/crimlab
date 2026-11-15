# Every plot is drawn on a null device; the data it returns are checked
# against the numbers of the result it drew.
drw <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expr
}

# the panel data that a time-adoption design draws
did <- function() {
  d <- expand.grid(year = 2015:2024, town = LETTERS[1:12])
  i <- match(d$town, LETTERS)
  d$first <- rep(c(2019, 2021, NA), each = 4)[i]
  d$bylaw <- as.integer(!is.na(d$first) & d$year >= d$first)
  d$rate <- 30 + i + 0.5 * (d$year - 2015) - 4 * d$bylaw + sin(seq_len(120))
  d
}

# the number of pages an expression draws
pages <- function(expr) {
  f <- tempfile("pg")
  grDevices::pdf(paste0(f, "%03d.pdf"), onefile = FALSE)
  on.exit(unlink(Sys.glob(paste0(f, "*"))), add = TRUE)
  expr
  grDevices::dev.off()
  length(Sys.glob(paste0(f, "*.pdf")))
}
