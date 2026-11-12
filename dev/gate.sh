#!/bin/sh
# The local gate: what CI checks, in one command, from the package root.
#   sh dev/gate.sh [LIB]        LIB: a private library to install into
#   GATE="test lint" sh dev/gate.sh    run only some steps
# Steps: doc test lint cov gp check lean kern. Each prints one GATE line.
# kern: every kernel family with a standalone test (src/t_<f>.cpp using
# tk.h, next to src/<f>.cpp) built and run without R at -O2, -Os,
# -O2 -ffp-contract=off and under ASan/UBSan; KBENCH=1 adds the timings.
set -u
lib=${1:-$HOME/tmp/cllib}
steps=${GATE:-doc test lint cov check}
gp=${GPLIB:-$HOME/tmp/gplib}
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
    Rscript -e 'roxygen2::roxygenise()' >/dev/null 2>&1 && say doc ok ||
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
    say cov "$out" ;;
  gp)
    out=$(Rscript -e '.libPaths(c(commandArgs(TRUE)[1], .libPaths()))
      ck <- goodpractice::all_checks(); ck <- ck[!grepl("^(covr|rcmdcheck)", ck)]
      f <- goodpractice::failed_checks(goodpractice::gp(".", checks = ck, quiet = TRUE))
      cat(length(ck), "checks, failed:", length(f), f)' "$gp" 2>/dev/null | tail -1)
    clean
    say gp "$out" ;;
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
  kern)
    t=$(mktemp -d); nb=0; bad=""
    for ts in src/t_*.cpp; do
      f=${ts#src/t_}; f=${f%.cpp}
      [ -f "src/$f.cpp" ] && grep -q '"tk.h"' "$ts" || continue
      for fl in "-O2" "-Os" "-O2 -ffp-contract=off" \
        "-O1 -g -fsanitize=address,undefined -fno-sanitize-recover=all"; do
        nb=$((nb + 1))
        # shellcheck disable=SC2086
        if ${CXX:-g++} -std=c++17 -Wall -Wextra -pedantic -Werror \
          -fno-exceptions -fno-rtti -DCLAB_MAIN $fl "src/$f.cpp" "$ts" \
          src/t_hmain.cpp -o "$t/th" >"$t/log" 2>&1 && "$t/th" >>"$t/log" 2>&1
        then
          tail -1 "$t/log" | sed "s/^/  $f [$fl] /"
        else
          cat "$t/log"; bad="$bad $f[$fl]"
        fi
      done
      [ "${KBENCH:-0}" = 1 ] && ${CXX:-g++} -std=c++17 -O2 -DCLAB_MAIN \
        "src/$f.cpp" "$ts" src/t_hmain.cpp -o "$t/tb" && nice -n 10 "$t/tb" bench
    done
    rm -rf "$t"
    say kern "$nb builds, failed:${bad:- none}"
    [ -z "$bad" ] || rc=1 ;;
  *) say "$s" "unknown step"; rc=1 ;;
  esac
done
exit $rc
