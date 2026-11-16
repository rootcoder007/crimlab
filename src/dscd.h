// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The text work of data description between the arguments and the
// results, shared by the R binding (bdscd.cpp) and clab_jrun()
// (cli/jdsc.inc): the cleaned column names, the missingness map, the
// missingness patterns and the selection of the groups a distinctions test
// compares. Nothing here knows about R or JSON. Text is UTF-8; labels are
// sorted in byte order. The numeric tables are clab_spark(),
// clab_prof(), clab_mhtab(), clab_mruns() and clab_topc() (dscd.cpp).
#ifndef CLAB_DSCD_H
#define CLAB_DSCD_H

#include <cstddef>
#include <string>
#include <vector>

namespace dscd {

// clncolnms(): the names nm folded to ASCII (clab_ascii()), "%", "&" and
// "#" spelled out, split into words at every run of other characters and
// before a capital that follows a small letter or a digit, joined in case
// cs (0 snake, 1 lower camel, 2 upper camel, 3 screaming snake, 4 as they
// are) with sep, "x" before a leading digit and for an empty name; a name
// repeated is numbered from its second time (sep 2, sep 3, ...).
std::vector<std::string> cnames(const std::vector<std::string> &nm, int cs, const std::string &sep);

// missgrid(): the lines of the map of na (n x p, nonzero missing) whose
// columns are named nm: a rule, the names cut to width characters written
// down, at most height row blocks shaded in quarters, the legend and a
// rule; Unicode blocks when uni, ASCII marks otherwise.
std::vector<std::string> mgrid(const int *na, std::size_t n, std::size_t p, const std::vector<std::string> &nm,
                               int height, int width, bool uni);

// misspat(): the missingness patterns of na (n x p), one 'X' (missing) or
// '.' per column, the maxpat most frequent (ties in pattern order), their
// rows, percent of rows, missing columns and the names of those columns.
struct pat {
  std::vector<std::string> code, cols;
  std::vector<int> rows, nmiss;
  std::vector<double> pct;
};
pat mpat(const int *na, std::size_t n, std::size_t p, const std::vector<std::string> &nm, int maxpat);

// dstntst(): of n groups g (gna: missing) with estimates y and interval
// lo, hi at the normal quantile z, the standard errors s = (hi - lo) / 2z,
// the estimates outside their own interval set to NaN (their groups in
// bad, once each, bna marking a missing one), the quality flags q (trimmed;
// none, or missing and empty, are "ok") not among ok, E, F, x, .., ... in
// odd and taken as ok, and use: the groups compared, neither suppressed
// (F, x, .., ...) nor pooling others (those named in pooled when hasp,
// else those whose label starts with total, indigenous identity,
// aboriginal identity or indigenous population, in any case) and with a
// finite estimate and a positive standard error.
struct sel {
  std::vector<char> use, bna;
  std::vector<double> y, s;
  std::vector<std::string> bad, odd;
};
sel dtsel(const std::vector<std::string> &g, const std::vector<char> &gna, const double *y, const double *lo,
          const double *hi, const std::vector<std::string> &q, const std::vector<char> &qna, bool hasp,
          const std::vector<std::string> &pooled, double z);

// mcartst(): the note on a test with nothing to test (df, o[1], at most
// 0) or with patterns skipped (o[7] of them), of clab_mcar()'s out; empty
// when there is neither.
std::string mcnote(const double *o);

}  // namespace dscd

#endif
