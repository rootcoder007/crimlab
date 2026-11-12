# Reference values and defect reproducer for src/aor.cpp (clab_aor).
# Field package: aoristic 2.0.0 (CRAN), R 4.5.0:
#   R_LIBS=/mnt/nvme/work/csls/lib R_ENVIRON_USER=/dev/null \
#   R_PROFILE_USER=/dev/null nice -n 10 /usr/bin/Rscript dev/doss/aoristic.R
suppressMessages(library(aoristic))
cat("aoristic", format(packageVersion("aoristic")), "\n")
f17 <- function(x) paste(sprintf("%.17g", x), collapse = ", ")
q <- function(e) suppressMessages(e)

# deterministic sample, rebuilt in src/t_aor.cpp (i = 0, ..., 299), UTC
i <- 0:299
from <- 1672617600 + 60 * ((i * 7919) %% 200000)
span <- 2 + (i * 104729) %% 480                       # minutes
span[i %% 50 == 7] <- 8 * 1440                        # longer than a week
span[i %% 50 == 21] <- -30                            # reversed
span[i %% 50 == 33] <- 1                              # one minute
to <- from + 60 * span
to[i %% 50 == 13] <- NA                               # no end
d <- data.frame(x = 0, y = 0,
  from = as.POSIXct(from, origin = "1970-01-01", tz = "UTC"),
  to = as.POSIXct(to, origin = "1970-01-01", tz = "UTC"))
r <- q(aoristic.df(d, "x", "y", "from", "to"))
h <- as.matrix(r[grep("^hour", names(r))])
cat("\n// colSums of aoristic.df(d, \"x\", \"y\", \"from\", \"to\") hour columns\n")
cat(sprintf("const double ref_tot[168] = {%s};\n", f17(colSums(h))))
cat("// rows 1:3 (i = 0, 1, 2)\n")
for (k in 1:3) cat(sprintf("const double ref_row%d[168] = {%s};\n", k - 1, f17(h[k, ])))

# defect (dossier aoristic.md section 3): daylight saving time
dst <- function(a, b) {
  e <- data.frame(x = 0, y = 0, from = as.POSIXct(a, tz = "Europe/London"),
                  to = as.POSIXct(b, tz = "Europe/London"))
  o <- q(aoristic.df(e, "x", "y", "from", "to"))
  w <- unlist(o[grep("^hour", names(o))])
  cat(a, "to", b, "(Europe/London): nonzero", names(w)[w > 0], w[w > 0],
      "| from, to as UTC seconds:", as.numeric(as.POSIXct(a, tz = "Europe/London")),
      as.numeric(as.POSIXct(b, tz = "Europe/London")), "\n")
}
cat("\n# spring forward (hour1 = Sunday 00h); truth Sunday 00h .25, 02h .5, 03h .25\n")
dst("2023-03-26 00:30", "2023-03-26 03:30")
cat("# fall back; truth Sunday 00h .125, 01h .5, 02h .25, 03h .125\n")
dst("2023-10-29 00:30", "2023-10-29 03:30")
cat("# Europe/London 2023 transitions (UTC seconds):",
    as.numeric(as.POSIXct("2023-03-26 01:00", tz = "UTC")),
    as.numeric(as.POSIXct("2023-10-29 01:00", tz = "UTC")), "\n")
cat("\n# n = 0:\n"); print(tryCatch(q(aoristic.df(d[0, ], "x", "y", "from", "to")), error = conditionMessage))

# timing: the dossier's case, n = 10,000 with spans up to 8 h, 1 run
set.seed(3); n <- 1e4
f <- as.POSIXct("2023-01-02 00:00", tz = "UTC") + runif(n, 0, 3e7)
big <- data.frame(x = runif(n), y = runif(n), from = f, to = f + runif(n, 0, 8 * 3600))
cat(sprintf("\nTIME aoristic.df n=10000 (8 h spans): %.3f s (1 run)\n",
  system.time(q(aoristic.df(big, "x", "y", "from", "to")))[["elapsed"]]))
