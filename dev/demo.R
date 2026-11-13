# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
# The one-minute demo behind man/figures/demo.svg: a veil-of-darkness table,
# a reported odds ratio checked against every labelling of its table, a
# clearance paradox, a categorical audit, a posterior from a likelihood
# ratio, hot spots, spatial autocorrelation, the command-line program, and
# the theorems behind a function. Every command really runs in this R
# session; only the prompt and the typing are paced for reading.
#
#   Rscript dev/demo.R              play it in the terminal
#   Rscript dev/demo.R demo.cast    write an asciicast (version 2) instead,
#                                   with no terminal recorder needed
# and then, for the picture in the README:
#   npx svg-term-cli@2.1.1 --in demo.cast --out man/figures/demo.svg \
#     --window --no-cursor --padding 12
suppressPackageStartupMessages(library(crimlab))
options(width = 88)

a <- commandArgs(TRUE)
cast <- if (length(a)) a[1L] else ""
wide <- 92L
high <- 27L
now <- 0
ev <- character(0)

esc <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\033", "\\u001b", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  gsub("\n", "\\n", s, fixed = TRUE)
}
# show text, then wait dt seconds
out <- function(s, dt = 0) {
  if (nzchar(cast)) {
    ev[length(ev) + 1L] <<- sprintf("[%.3f, \"o\", \"%s\"]", now,
                                    esc(gsub("\n", "\r\n", s, fixed = TRUE)))
    now <<- now + dt
  } else {
    cat(s)
    Sys.sleep(dt)
  }
}
say <- function(s) out(sprintf("\033[2m# %s\033[0m\n", s), 0.7)
# type a shell command at a $ prompt and run it (the command-line program is
# cli/crimlab, or the path in CLAB_CLI)
sh <- function(cmd, pause = 1.4) {
  out("\033[1;32m$\033[0m ")
  ch <- strsplit(cmd, "", fixed = TRUE)[[1L]]
  for (i in seq(1L, length(ch), by = 3L)) {
    out(paste(ch[i:min(i + 2L, length(ch))], collapse = ""), 0.07)
  }
  out("\n", 0.35)
  exe <- Sys.getenv("CLAB_CLI", "cli/crimlab")
  res <- system(sub("^crimlab", shQuote(exe), cmd), intern = TRUE)
  out(paste0(paste(res, collapse = "\n"), "\n"), pause)
}
# type a command at the prompt, run it, print what R would print
run <- function(code, pause = 1.4) {
  out("\033[1;32m>\033[0m ")
  ch <- strsplit(code, "", fixed = TRUE)[[1L]]
  for (i in seq(1L, length(ch), by = 3L)) {
    out(paste(ch[i:min(i + 2L, length(ch))], collapse = ""), 0.07)
  }
  out("\n", 0.35)
  v <- withVisible(eval(parse(text = code), globalenv()))
  if (v$visible) {
    out(paste0(paste(utils::capture.output(print(v$value)), collapse = "\n"), "\n"))
  }
  out("", pause)
}

say("stops by group, in daylight and after dark")
run("day <- c(Black = 130, White = 220); night <- c(Black = 95, White = 230)", 0.4)
say("veil of darkness: each group's share of stops, lit and unlit")
run("vodtst(day, night)", 2.2)
say("a report gives an odds ratio of 0.44 for A against B. Which way round is it?")
run("tab <- rbind(A = c(yes = 20, no = 80), B = c(yes = 10, no = 90))", 0.4)
run("relchk(tab, 0.44, ref = \"B\", event = \"yes\")", 3)
say("clearance fell for violent crime and for property crime; the total rose")
run("a <- list(clr = c(violent = 60, property = 180), rec = c(100, 900))", 0.3)
run("b <- list(clr = c(violent = 150, property = 105), rec = c(300, 700))", 0.3)
run("round(clrsmp(a$clr, a$rec, b$clr, b$rec)$total, 3)", 3)
say("what would an analysis make of this column?")
run("catchk(data.frame(sex = c(\"Male\", \"male\", \"Female\", \"Female\", NA)), nmin = 2)", 3)
say("a likelihood ratio of 1000 on prior odds of 1 in 10,000: the posterior probability")
run("lrpost(1000, prior = 1 / 10000)", 2)
say("hot spots: nearest-neighbour hierarchical clusters of 120 incidents")
run("set.seed(2); x <- c(rnorm(40, 0, .3), rnorm(40, 6, .3), runif(40, -5, 10))", 0.3)
run("y <- c(rnorm(40, 0, .3), rnorm(40, 6, .3), runif(40, -5, 10))", 0.3)
run("round(nnhclu(x, y, area = 225, minpts = 5)$clusters[, c(\"points\", \"mx\", \"my\")], 2)", 2.5)
say("are neighbouring places alike? Geary's C on a 4-nearest-neighbour graph")
run("round(gearc(x + y, spwt(cbind(x, y), k = 4)), 3)", 2.5)
say("the same kernels from the shell")
sh("crimlab lrpost --lr 1000 --prior 0.0001", 2)
say("the machine-checked theorems behind a function")
run("proof(\"vodtst\")", 2.5)
out(" ", 5)

if (nzchar(cast)) {
  writeLines(c(sprintf("{\"version\": 2, \"width\": %d, \"height\": %d}", wide, high),
               ev, sprintf("[%.3f, \"o\", \" \"]", now)), cast)
  cat(length(ev), "events,", round(now, 1), "seconds ->", cast, "\n")
}
