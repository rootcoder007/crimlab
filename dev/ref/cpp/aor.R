# Reference values for src/aor.cpp (clab_aor).
# Field package: aoristic 2.0.0 (CRAN), R 4.5.0.
source("dev/ref/cpp/pre.R")
suppressMessages(library(aoristic))
cat("aoristic", format(packageVersion("aoristic")), "\n")
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

# the Europe/London 2023 clock changes used by the test (UTC seconds)
cat("// transitions:", as.numeric(as.POSIXct("2023-03-26 01:00", tz = "UTC")),
    as.numeric(as.POSIXct("2023-10-29 01:00", tz = "UTC")), "\n")

# timing: n = 10,000 with spans up to 8 h, 1 run
set.seed(3); n <- 1e4
f <- as.POSIXct("2023-01-02 00:00", tz = "UTC") + runif(n, 0, 3e7)
big <- data.frame(x = runif(n), y = runif(n), from = f, to = f + runif(n, 0, 8 * 3600))
cat(sprintf("\nTIME aoristic.df n=10000 (8 h spans): %.3f s (1 run)\n",
  system.time(q(aoristic.df(big, "x", "y", "from", "to")))[["elapsed"]]))
