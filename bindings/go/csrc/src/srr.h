// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// rOpenSci statistical software standards (General and EDA) that hold for
// the package as a whole, read by the srr package (srr::srr_report()).
// Standards met in one place are tagged there: the checks in rchk.cpp and
// bind.cpp, the kernels in core.cpp, the self-test in selftest.cpp and
// t_*.cpp, the R tests, the plots' help pages and the vignette. srr reads
// //' comments in src/; a block is read only when a line of code follows
// it. Nothing here is compiled: no source file includes this header.
#ifndef CLAB_SRR_H
#define CLAB_SRR_H

//' srr_stats
//'
//' Standards the package meets.
//'
//' @srrstats {G1.1} README "Algorithms" says that every method is a
//'   published one, none new, and that crimlab claims no first
//'   implementation in R.
//' @srrstats {G1.2} README "Life cycle" states the stage: experimental,
//'   so that names, arguments and results may change before version 1.0,
//'   with each change recorded in NEWS.md.
//' @srrstats {G1.4a} Every internal R function (the ones named with a
//'   leading dot) has a roxygen2 comment saying what it does, ending in
//'   `@noRd`.
//' @srrstats {G5.10} The whole suite is the extended set: tests/testthat.R
//'   runs every test file when the environment variable NOT_CRAN is
//'   "true" (the project's workflows set it, as devtools does) and one
//'   file otherwise, which is what CRAN runs.
//' @srrstats {G5.12} tests/README.md gives the switch, what the full
//'   suite needs and how long it takes.
//' @srrstats {G1.6} The package makes no claim of speed or accuracy over
//'   other packages; tests named "agree with" (test-disp.R, test-pol.R,
//'   test-corr.R) compare its results with another R implementation when
//'   that package is installed. The speed of enpath(), rfor()
//'   and dmlx() against third-party packages is measured, not claimed:
//'   dev/bench-dml.R checks that the estimates agree, then times both
//'   sides and writes inst/bench/dml.csv with the machine, its load and
//'   the versions.
//' @srrstats {G2.0a} Each @param states the length expected: "one number",
//'   "same length", "one per locus", "once or once per ratio".
//' @srrstats {G2.1a} Each @param states the type expected: numbers,
//'   TRUE or FALSE, allele labels, names of `p`, function names.
//' @srrstats {G2.14} Missing data are always an error (G2.14a); no option
//'   drops them, for the reasons under G2.14b.
//' @srrstats {G2.14a} An NA in any argument stops with an error naming it,
//'   except where a help page says the value is not read (sellab()'s y for
//'   detained defendants, gtmp()'s q at a homozygous locus).
//' @srrstats {G5.2} The refusals of every exported function are triggered
//'   by the errs() tables of test-disp.R, test-frns.R, test-frns2.R,
//'   test-pol.R and test-corr.R, and by test-chk.R and test-bind.R; the
//'   kernels' own error codes by the self-test (guards in src/t_*.cpp).
//' @srrstats {G5.4b} Tests named "agree with" compare results with an
//'   earlier R implementation of the same methods when it is installed.
//' @srrstats {G5.8} Edge conditions are tested; see G5.8a-d.
//' @srrstats {G5.8a} Zero-length input is refused (losbnd() and lrcal() in
//'   test-chk.R, test-corr.R and test-frns2.R; the kernels in the guards of
//'   src/t_*.cpp).
//' @srrstats {G5.8b} Character input to numeric arguments is refused
//'   (losbnd("1") in test-chk.R, outtst("a") in test-disp.R).
//' @srrstats {G5.8d} Input outside a method's scope is refused: chances
//'   above 1 (test-disp.R), frequencies that do not form a distribution
//'   (test-edge.R, test-frns2.R), more than 20 unexplained alleles in a
//'   mixture (src/t_frns.cpp).
//' @srrstats {G5.9} Noise susceptibility is tested under G5.9a
//'   (src/t_disp.cpp) and G5.9b (ntlvn()'s seed, test-net.R).
//' @srrstats {EA1.1} README "What it takes" lists the kinds of data: counts
//'   by group, individual binary outcomes, allele frequencies and
//'   genotypes, and likelihood ratios from validation studies.
//' @srrstats {EA1.3} Each @param gives the kind of data the argument takes.
//' @srrstats {EA4.2} Results are data frames, named vectors and lists of
//'   these. The lists that are fits or reports carry a class with a
//'   print method of its own; those that are drawn have a plot method
//'   (the help page `plotm` lists most of them, the pages of their
//'   functions the rest), and some have a summary method.
//' @srrstats {EA6.0} Return values of every exported function are tested;
//'   see EA6.0a-e.
//' @srrstats {EA6.0e} Nearly every expect_equal() on a number, and every
//'   tk_eq() and tk_veq() of the self-test, gives an explicit tolerance;
//'   the few R ones that do not use testthat edition 3's default, the
//'   square root of machine precision.
//' @srrstats {EA6.1} The data each plot draws is returned and tested
//'   ("the plots return what they draw" in test-frns2.R).
//' @srrstats {EA3.0} enpath() chooses the penalty of the elastic net by
//'   cross-validation over the caller's folds and returns the whole path
//'   with the losses, so the choice can be checked.
//' @srrstats {G5.9b} rfor() draws its bags and columns from a seed:
//'   test-rf.R checks that one seed gives one forest, that another gives
//'   another, and that the error stays within its bound across seeds.
//' @noRd
#define CLAB_SRR_MET 1

//' NA_standards
//'
//' Standards that do not apply, with the reason for each.
//'
//' @srrstatsNA {G1.5} There is no associated publication with performance
//'   claims to reproduce.
//' @srrstatsNA {G2.4d} No argument has to be a factor, so no input needs
//'   converting to one; where a function counts over a fixed set of
//'   levels it builds the factor itself, inside table().
//' @srrstatsNA {G2.14b} Dropping a missing count would change the groups or
//'   the population a result describes, so missing values are an error,
//'   never silently ignored; the caller can subset first.
//' @srrstatsNA {G2.14c} Imputing a count, a frequency or an outcome would
//'   invent the data the result is about.
//' @srrstatsNA {G5.0} No standard data set (such as NIST's) exists for these
//'   quantities; tests use small inputs with answers computed by hand.
//' @srrstatsNA {G5.1} The package ships no data set; test data are built
//'   in the test files, and the files tests/testthat/*.rds hold reference
//'   values that scripts in dev/ write from other packages (netref.R,
//'   ref-dml.R, ref-fair.R, ref-learn.R, ref-rf.R).
//' @srrstatsNA {G5.11} No test needs a large data set or a download: the
//'   tests of the open tables read a copy of the service that
//'   helper-od.R makes in a temporary folder.
//' @srrstatsNA {G5.11a} No test downloads anything (G5.11).
//' @srrstatsNA {EA2.0} No table filtering or joining.
//' @srrstatsNA {EA2.1} No index columns, because no joins (EA2.0).
//' @srrstatsNA {EA2.2} No index columns (EA2.0).
//' @srrstatsNA {EA2.2a} No index columns (EA2.0).
//' @srrstatsNA {EA2.2b} No index columns (EA2.0).
//' @srrstatsNA {EA2.3} No table joins.
//' @srrstatsNA {EA2.4} No multi-table input.
//' @srrstatsNA {EA2.5} No index columns (EA2.0).
//' @srrstatsNA {EA3.1} Each function computes one named statistic; there are
//'   no alternative processes or models to compare in a standard form.
//' @srrstatsNA {EA5.1} No plot sets a typeface.
//' @srrstatsNA {EA5.6} No visualisation library is bundled; the plots use
//'   base graphics.
//' @noRd
#endif
