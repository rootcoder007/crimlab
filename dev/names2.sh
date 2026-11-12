#!/bin/sh
# Name gate for the kernels: every clab_<name> declared in src/crimlab.h
# becomes the R function <name>(), so <name> must be free.
#   BASE=main.h sh dev/names2.sh [HEADER] [OTHER_HEADER ...]
# (main.h: e.g. git show origin/main:src/crimlab.h > main.h)
# Checks, for each bare name: (i) no CRAN package exports a function of
# that name (pkgcheck's table, the one the rOpenSci bot uses); (ii) it is
# not defined with base R and the recommended packages stats, utils,
# graphics, methods, MASS, survival, spatial, cluster, nlme attached;
# (iii) no CRAN package has that name; (iv) not exported by rmoriebricklayer,
# rmoriedata or rmorie (NAMESPACE files under $WORK, default ~/work);
# (v) not declared in another branch's header unless both inherit it from
# the base header BASE. Needs pkgcheck installed and network
# access for the CRAN package list. Prints one line per clash and exits 1
# if there is any.
set -u
hdr=${1:-src/crimlab.h}
[ $# -gt 0 ] && shift
work=${WORK:-$HOME/work}
R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null
export R_ENVIRON_USER R_PROFILE_USER
names=$(grep -h -o 'clab_[a-z0-9]*(' "$hdr" | sed 's/^clab_//; s/($//' |
  grep -v '^strerror$' | sort -u)
base=${BASE:-/dev/null}
dups=$( (grep -h -o 'clab_[a-z0-9]*(' "$hdr" | sort -u
  for o in "$@"; do grep -h -o 'clab_[a-z0-9]*(' "$o" | sort -u; done) |
  sort | uniq -d |
  while read -r d; do grep -q "$d" "$base" || echo "$d"; done |
  sed 's/^clab_//; s/($//')
Rscript -e '
a <- commandArgs(TRUE); work <- a[1]; dups <- strsplit(a[2], " ")[[1]]
x <- a[-(1:2)]
bad <- character()
f <- get("fn_names_on_cran", asNamespace("pkgcheck"))(x)
if (nrow(f)) bad <- c(bad, sprintf("%s: CRAN function in %s", f$fn_name, f$package))
for (p in c("stats", "utils", "graphics", "methods", "MASS", "survival",
  "spatial", "cluster", "nlme")) suppressMessages(library(p, character.only = TRUE))
e <- x[vapply(x, exists, NA)]
if (length(e)) bad <- c(bad, sprintf("%s: defined in base or recommended R", e))
cp <- rownames(available.packages(repos = "https://cloud.r-project.org"))
k <- x[tolower(x) %in% tolower(cp)]
if (length(k)) bad <- c(bad, sprintf("%s: a CRAN package name", k))
for (p in c("rmoriebricklayer", "rmoriedata", "rmorie")) {
  ns <- file.path(work, p, "NAMESPACE")
  if (!file.exists(ns)) { bad <- c(bad, sprintf("(%s: no NAMESPACE at %s)", p, ns)); next }
  ex <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", readLines(ns), value = TRUE))
  k <- intersect(x, gsub("\"", "", ex))
  if (length(k)) bad <- c(bad, sprintf("%s: exported by %s", k, p))
}
d <- dups[nzchar(dups)]
if (length(d)) bad <- c(bad, sprintf("%s: declared in two headers", d))
cat(sprintf("names checked: %d\n", length(x)))
if (length(bad)) { writeLines(bad); quit(status = 1) }
cat("all free\n")
' "$work" "$(echo $dups)" $names
