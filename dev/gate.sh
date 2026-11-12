#!/bin/sh
# The local gate: what CI checks, in one command, from the package root.
#   sh dev/gate.sh [LIB]        LIB: a private library to install into
#   GATE="test lint" sh dev/gate.sh    run only some steps
# Steps: doc test lint cov gp check lean srr vig site. Each prints one
# GATE line and sets a non-zero exit on failure (cov: below COVMIN, 100 by
# default; gp: any failed check but urlchecker_ok). doc, gp and srr need the
# srr package (the roxygen roclet in DESCRIPTION): SRRLIB names a private
# library holding it.
set -u
lib=${1:-$HOME/tmp/cllib}
steps=${GATE:-doc test lint cov check}
gp=${GPLIB:-$HOME/tmp/gplib}
srrlib=${SRRLIB:-$HOME/tmp/srrlib}
lean=${LEANPK:-$HOME/work/rmorie/research/lean/.lake/packages}
R_ENVIRON_USER=/dev/null
R_PROFILE_USER=/dev/null
export R_ENVIRON_USER R_PROFILE_USER
mkdir -p "$lib"
rc=0
say() { printf 'GATE %s: %s\n' "$1" "$2"; }
clean() { rm -f src/*.o src/*.so; }
for s in $steps; do
  case $s in
  doc)
    Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      roxygen2::roxygenise()' "$srrlib" >/dev/null 2>&1 && say doc ok ||
      { say doc FAILED; rc=1; }
    clean ;;
  test)
    R CMD INSTALL --library="$lib" . >/dev/null 2>&1 || { say install FAILED; rc=1; }
    clean
    out=$(Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      r <- as.data.frame(testthat::test_dir("tests/testthat", package = "crimlab",
        load_package = "installed", reporter = "silent", stop_on_failure = FALSE))
      cat(sum(r$passed), sum(r$failed), sum(r$skipped), sum(r$error))' "$lib" 2>/dev/null | tail -1)
    say test "pass fail skip error = $out"
    case $out in *" 0 0 0") ;; *) rc=1 ;; esac ;;
  lint)
    n=$(Rscript -e 'cat(length(lintr::lint_package()))' 2>/dev/null | tail -1)
    say lint "$n lints"; [ "$n" = 0 ] || rc=1 ;;
  cov)
    out=$(Rscript -e 'x <- covr::package_coverage(".")
      z <- covr::tally_coverage(x); z <- z[z$value == 0, ]
      cat(sprintf("%.2f%% uncovered: %s", covr::percent_coverage(x),
        paste(sprintf("%s:%d", z$filename, z$line), collapse = " ")))' 2>/dev/null | tail -1)
    clean
    say cov "$out"
    # below COVMIN per cent (default 100), or no number at all, fails
    pct=${out%%%*}
    awk -v p="$pct" -v m="${COVMIN:-100}" 'BEGIN { exit !(p + 0 >= m + 0 && p != "") }' ||
      rc=1 ;;
  gp)
    # srr defines the @srrstats tags; without it roxygen calls them unknown.
    # Any failed check but urlchecker_ok (a site not yet published) fails.
    out=$(Rscript -e '.libPaths(c(commandArgs(TRUE)[1:2], .libPaths()))
      loadNamespace("srr")
      ck <- goodpractice::all_checks(); ck <- ck[!grepl("^(covr|rcmdcheck)", ck)]
      f <- goodpractice::failed_checks(goodpractice::gp(".", checks = ck, quiet = TRUE))
      cat(length(ck), "checks, failed:", length(f), f)' "$gp" "$srrlib" 2>/dev/null | tail -1)
    clean
    say gp "$out"
    case $out in *"failed: 0"|*"failed: 1 urlchecker_ok") ;; *) rc=1 ;; esac ;;
  check)
    t=$(mktemp -d)
    (cd "$t" && R CMD build "$OLDPWD" >/dev/null 2>&1 &&
      _R_CHECK_FORCE_SUGGESTS_=false R CMD check --as-cran --no-manual \
        crimlab_*.tar.gz 2>&1 | grep '^Status' | tail -1) > "$t/st" 2>/dev/null
    st=$(cat "$t/st"); rm -rf "$t"
    say check "${st:-FAILED}"
    case $st in "Status: OK"|"Status: 1 NOTE"|"Status: 2 NOTEs") ;; *) rc=1 ;; esac ;;
  lean)
    t=$(mktemp -d); mkdir "$t/.lake"; cp inst/lean/* "$t/"
    ln -s "$lean" "$t/.lake/packages"
    st=$(cd "$t" && "$HOME/.elan/bin/lake" build 2>&1 | tail -1)
    so=$(grep -l sorry inst/lean/*.lean 2>/dev/null | wc -l)
    rm -rf "$t"
    say lean "$st; files with sorry: $so"
    case $st in "Build completed successfully"*) [ "$so" = 0 ] || rc=1 ;; *) rc=1 ;; esac ;;
  srr)
    # every standard of the declared categories tagged, none left TODO
    out=$(Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      m <- srr::srr_report(".", view = FALSE, roxygenise = FALSE)
      n <- function(tag) {
        i <- grep(paste0("^## Standards with `", tag, "` tag"), m)
        if (!length(i)) return("0")
        x <- m[seq(i, length(m))]
        sub(".*: ([0-9]+) / .*", "\\1", grep("^- +Total", x, value = TRUE)[1])
      }
      bad <- any(grepl("heavy_multiplication_x|^## Missing Standards", m))
      cat(sprintf("complied %s NA %s TODO %s%s", n("srrstats"),
        n("srrstatsNA"), n("srrstatsTODO"), if (bad) " ERRORS" else ""))
      ' "$srrlib" 2>/dev/null | tail -1)
    say srr "$out"
    case $out in *ERRORS*|"") rc=1 ;; esac ;;
  vig)
    t=$(mktemp -d)
    Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      tools::buildVignettes(dir = ".", tangle = TRUE, quiet = TRUE)' "$lib" \
      >"$t/log" 2>&1 && say vig ok || { say vig FAILED; tail -5 "$t/log"; rc=1; }
    rm -rf "$t" vignettes/*.html vignettes/*.R vignettes/*.md; clean ;;
  site)
    # SITEDIR keeps the site (CI deploys it); otherwise it goes to a temp dir
    t=$(mktemp -d); dest=${SITEDIR:-$t/site}
    Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      pkgdown::build_site(".", override = list(destination = commandArgs(TRUE)[2]),
        preview = FALSE, install = FALSE, new_process = FALSE)' "$lib" "$dest" \
      >"$t/log" 2>&1
    w=$(grep -ci "warn\|missing" "$t/log")
    if [ -f "$dest/index.html" ] && ! grep -q "^Error" "$t/log"; then
      say site "built, $w lines mention a warning or missing topic"
      [ "$w" = 0 ] || { grep -i "warn\|missing" "$t/log"; rc=1; }
    else say site FAILED; tail -20 "$t/log"; rc=1; fi
    rm -rf "$t"; clean ;;
  *) say "$s" "unknown step"; rc=1 ;;
  esac
done
exit $rc
