// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of data description between arguments and results
// (src/dscd.cpp, dscd.h). The R wrappers (R/dsc.R) check the arguments and
// name and class the results; these take the columns apart, call the work
// that clab_jrun() calls too, and hand back plain vectors. A failed kernel
// throws, so no C++ object is alive when R's error is raised (glue.cpp).
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "bcol.h"
#include "dscd.h"
#include "rchk.h"

namespace {

using std::int64_t;
using std::size_t;

void ok(int rc) {
  clab::ichk(rc);
  if (rc != CLAB_OK) throw std::runtime_error(clab_strerror(rc));
}

// a character vector as UTF-8 strings, NA as "" flagged in na
std::vector<std::string> utf8(SEXP x, std::vector<char> *na = nullptr) {
  const R_xlen_t n = Rf_isNull(x) ? 0 : XLENGTH(x);
  std::vector<std::string> r(static_cast<size_t>(n));
  if (na != nullptr) na->assign(static_cast<size_t>(n), 0);
  for (R_xlen_t i = 0; i < n; ++i) {
    SEXP e = STRING_ELT(x, i);
    if (e == NA_STRING) {
      if (na != nullptr) (*na)[static_cast<size_t>(i)] = 1;
      continue;
    }
    r[static_cast<size_t>(i)] = Rf_translateCharUTF8(e);
  }
  return r;
}

// a character vector of UTF-8 strings s, those flagged in na NA
SEXP chr(const std::vector<std::string> &s, const std::vector<char> *na = nullptr) {
  SEXP v = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(s.size())));
  for (size_t i = 0; i < s.size(); ++i) {
    const bool miss = na != nullptr && (*na)[i];
    SET_STRING_ELT(v, static_cast<R_xlen_t>(i), miss ? NA_STRING : Rf_mkCharCE(s[i].c_str(), CE_UTF8));
  }
  UNPROTECT(1);
  return v;
}

SEXP ints(const std::vector<int> &x) {
  SEXP v = Rf_allocVector(INTSXP, static_cast<R_xlen_t>(x.size()));
  for (size_t i = 0; i < x.size(); ++i) INTEGER(v)[i] = x[i];
  return v;
}

SEXP dbls(const std::vector<double> &x) {
  SEXP v = Rf_allocVector(REALSXP, static_cast<R_xlen_t>(x.size()));
  for (size_t i = 0; i < x.size(); ++i) REAL(v)[i] = std::isnan(x[i]) ? NA_REAL : x[i];
  return v;
}

// a list of the SEXPs v, each protected by the caller
SEXP lst(std::initializer_list<SEXP> v) {
  SEXP r = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(v.size())));
  R_xlen_t i = 0;
  for (SEXP e : v) SET_VECTOR_ELT(r, i++, e);
  UNPROTECT(1);
  return r;
}

// the rows and columns of a logical matrix (a vector: one column)
void dims(SEXP na, size_t &n, size_t &p) {
  SEXP d = Rf_getAttrib(na, R_DimSymbol);
  const bool m = !Rf_isNull(d) && XLENGTH(d) == 2;
  n = m ? static_cast<size_t>(INTEGER(d)[0]) : static_cast<size_t>(XLENGTH(na));
  p = m ? static_cast<size_t>(INTEGER(d)[1]) : 1;
  clab::ilen(static_cast<double>(n), "data");
}

}  // namespace

extern "C" {

SEXP C_clab_r_spark(SEXP x, SEXP bins, SEXP uni) {
  clab::who = "inlhist";
  const int b = Rf_asInteger(bins);
  int64_t len = 0;
  int rc = CLAB_OK;
  std::string s;
  {
    std::vector<char> buf(static_cast<size_t>(b) * 3 + 1, 0);
    rc = clab_spark(REAL(x), clab::xlen(x, "x"), b, Rf_asLogical(uni) == TRUE, buf.data(), &len);
    s.assign(buf.data(), static_cast<size_t>(len));
  }
  ok(rc);
  return chr({s});
}

SEXP C_clab_r_cnames(SEXP nms, SEXP cs, SEXP sep) {
  clab::who = "clncolnms";
  // the bytes as they are, as toascii() reads them
  const R_xlen_t n = clab::xlen(nms, "data");
  std::vector<std::string> in(static_cast<size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) in[static_cast<size_t>(i)] = CHAR(STRING_ELT(nms, i));
  return chr(dscd::cnames(in, Rf_asInteger(cs), Rf_translateCharUTF8(STRING_ELT(sep, 0))));
}

SEXP C_clab_r_mgridl(SEXP na, SEXP nms, SEXP height, SEXP width, SEXP uni) {
  clab::who = "missgrid";
  size_t n = 0, p = 0;
  dims(na, n, p);
  return chr(dscd::mgrid(LOGICAL(na), n, p, utf8(nms), Rf_asInteger(height), Rf_asInteger(width),
                         Rf_asLogical(uni) == TRUE));
}

SEXP C_clab_r_mpat(SEXP na, SEXP nms, SEXP maxpat) {
  clab::who = "misspat";
  size_t n = 0, p = 0;
  dims(na, n, p);
  const dscd::pat r = dscd::mpat(LOGICAL(na), n, p, utf8(nms), Rf_asInteger(maxpat));
  SEXP a = PROTECT(chr(r.code));
  SEXP b = PROTECT(ints(r.rows));
  SEXP c = PROTECT(dbls(r.pct));
  SEXP d = PROTECT(ints(r.nmiss));
  SEXP e = PROTECT(chr(r.cols));
  SEXP o = lst({a, b, c, d, e});
  UNPROTECT(5);
  return o;
}

// profcols(): cols, a list of numeric vectors or NULL (not numeric); the
// statistics as columns and the p x np matrix of quantiles.
SEXP C_clab_r_prof(SEXP cols, SEXP qnts, SEXP uni) {
  clab::who = "profcols";
  static const char *const nm[] = {"n_distinct", "nzero", "nneg", "ninf", "mean", "sd", "median",
                                   "mad", "min", "max", "skewness", "n_outliers", "hist", "top"};
  static const bool isint[] = {true, true, true, true, false, false, false,
                               false, false, false, false, true};
  const R_xlen_t p = XLENGTH(cols);
  const int np = static_cast<int>(clab::xlen(qnts, "qnts"));
  SEXP st = PROTECT(Rf_allocVector(VECSXP, 14));
  for (int k = 0; k < 12; ++k) SET_VECTOR_ELT(st, k, Rf_allocVector(isint[k] ? INTSXP : REALSXP, p));
  SET_VECTOR_ELT(st, 12, Rf_allocVector(STRSXP, p));
  SET_VECTOR_ELT(st, 13, Rf_allocVector(STRSXP, p));
  SEXP qm = PROTECT(Rf_allocMatrix(REALSXP, static_cast<int>(p), np));
  const int u = Rf_asLogical(uni) == TRUE;
  for (R_xlen_t j = 0; j < p; ++j) {
    SEXP x = VECTOR_ELT(cols, j);
    double s[12] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int rc = CLAB_OK;
    std::string h;
    bool hist = false;
    {
      std::vector<double> q(static_cast<size_t>(np), 0.0);
      const int64_t n = Rf_isNull(x) ? 0 : clab::xlen(x, "data");
      rc = clab_prof(Rf_isNull(x) ? nullptr : REAL(x), n, REAL(qnts), np, s, q.data());
      for (int k = 0; k < np; ++k) {
        const double v = q[static_cast<size_t>(k)];
        REAL(qm)[static_cast<R_xlen_t>(k) * p + j] = std::isnan(v) ? NA_REAL : v;
      }
      if (rc == CLAB_OK && !std::isnan(s[8])) {
        std::vector<char> buf(31, 0);
        int64_t len = 0;
        rc = clab_spark(REAL(x), n, 10, u, buf.data(), &len);
        h.assign(buf.data(), static_cast<size_t>(len));
        hist = true;
      }
    }
    ok(rc);
    for (int k = 0; k < 12; ++k) {
      SEXP c = VECTOR_ELT(st, k);
      if (isint[k]) INTEGER(c)[j] = std::isnan(s[k]) ? NA_INTEGER : static_cast<int>(s[k]);
      else REAL(c)[j] = std::isnan(s[k]) ? NA_REAL : s[k];
    }
    SET_STRING_ELT(VECTOR_ELT(st, 12), j, hist ? Rf_mkCharCE(h.c_str(), CE_UTF8) : NA_STRING);
    SET_STRING_ELT(VECTOR_ELT(st, 13), j, NA_STRING);
  }
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 14));
  for (int k = 0; k < 14; ++k) SET_STRING_ELT(nms, k, Rf_mkChar(nm[k]));
  Rf_setAttrib(st, R_NamesSymbol, nms);
  SEXP o = lst({st, qm});
  UNPROTECT(3);
  return o;
}

// mahaoutl(): the table's columns, or 1 (too few complete rows) or 2
// (collinear columns).
SEXP C_clab_r_mhtab(SEXP x, SEXP robust, SEXP alpha) {
  clab::who = "mahaoutl";
  const int n = Rf_nrows(x), p = Rf_ncols(x);
  SEXP row = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP d = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP pv = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP lp = PROTECT(Rf_allocVector(REALSXP, n));
  SEXP ol = PROTECT(Rf_allocVector(LGLSXP, n));
  int why = 0;
  const int rc = clab_mhtab(REAL(x), n, p, Rf_asLogical(robust), Rf_asReal(alpha), INTEGER(row), REAL(d), REAL(pv),
                            REAL(lp), LOGICAL(ol), &why);
  if (why > 0) {
    UNPROTECT(5);
    return Rf_ScalarInteger(why);
  }
  ok(rc);
  for (int i = 0; i < n; ++i) {
    if (LOGICAL(ol)[i] < 0) {
      LOGICAL(ol)[i] = NA_LOGICAL;
      REAL(d)[i] = REAL(pv)[i] = REAL(lp)[i] = NA_REAL;
    }
  }
  SEXP o = lst({row, d, pv, lp, ol});
  UNPROTECT(5);
  return o;
}

// missruns(): the column (from 1), first and last row of every run.
SEXP C_clab_r_mruns(SEXP na, SEXP minrun, SEXP rank) {
  clab::who = "missruns";
  size_t n = 0, p = 0;
  dims(na, n, p);
  const size_t cap = p * (n / 2 + 1);
  std::vector<int> c(cap, 0), s(cap, 0), e(cap, 0);
  int64_t k = 0;
  const int rc = clab_mruns(LOGICAL(na), static_cast<int64_t>(n), static_cast<int>(p), Rf_asInteger(minrun),
                            INTEGER(rank), c.data(), s.data(), e.data(), &k);
  ok(rc);
  c.resize(static_cast<size_t>(k));
  s.resize(static_cast<size_t>(k));
  e.resize(static_cast<size_t>(k));
  SEXP a = PROTECT(ints(c));
  SEXP b = PROTECT(ints(s));
  SEXP f = PROTECT(ints(e));
  SEXP o = lst({a, b, f});
  UNPROTECT(3);
  return o;
}

// topcors(): the pairs (columns from 1) and their correlations, at most top
// of them, the strongest first.
SEXP C_clab_r_topcor(SEXP x, SEXP pear, SEXP minabs, SEXP top) {
  clab::who = "topcors";
  const int n = Rf_nrows(x), p = Rf_ncols(x);
  const size_t k = static_cast<size_t>(p) * static_cast<size_t>(p - 1) / 2;
  std::vector<double> r(k, 0.0);
  std::vector<int> np(k, 0), idx(k, 0), a, b;
  std::vector<double> v;
  clab_pcor(REAL(x), n, p, Rf_asLogical(pear), 3, r.data(), np.data());
  for (size_t i = 0; i < k; ++i) {
    if (np[i] < 3) r[i] = NAN;
  }
  int64_t m = 0;
  const int rc = clab_topc(r.data(), static_cast<int64_t>(k), Rf_asReal(minabs),
                           static_cast<int64_t>(Rf_asInteger(top)), idx.data(), &m);
  ok(rc);
  std::vector<int> pi, pj;
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      pi.push_back(i + 1);
      pj.push_back(j + 1);
    }
  }
  for (int64_t t = 0; t < m; ++t) {
    const size_t q = static_cast<size_t>(idx[static_cast<size_t>(t)] - 1);
    a.push_back(pi[q]);
    b.push_back(pj[q]);
    v.push_back(r[q]);
  }
  SEXP sa = PROTECT(ints(a));
  SEXP sb = PROTECT(ints(b));
  SEXP sv = PROTECT(dbls(v));
  SEXP o = lst({sa, sb, sv});
  UNPROTECT(3);
  return o;
}

// dstntst(): the groups compared (use), the estimates (those outside their
// interval NA), the standard errors, and the groups and quality flags to
// warn of.
SEXP C_clab_r_dtsel(SEXP g, SEXP y, SEXP lo, SEXP hi, SEXP q, SEXP pooled, SEXP z) {
  clab::who = "dstntst";
  std::vector<char> gna, qna;
  const std::vector<std::string> gs = utf8(g, &gna), qs = utf8(q, &qna), ps = utf8(pooled);
  const dscd::sel r = dscd::dtsel(gs, gna, REAL(y), REAL(lo), REAL(hi), qs, qna, !Rf_isNull(pooled), ps,
                                  Rf_asReal(z));
  SEXP use = PROTECT(Rf_allocVector(LGLSXP, static_cast<R_xlen_t>(r.use.size())));
  for (size_t i = 0; i < r.use.size(); ++i) LOGICAL(use)[i] = r.use[i];
  SEXP sy = PROTECT(dbls(r.y));
  SEXP ss = PROTECT(dbls(r.s));
  SEXP sb = PROTECT(chr(r.bad, &r.bna));
  SEXP so = PROTECT(chr(r.odd));
  SEXP o = PROTECT(lst({use, sy, ss, sb, so}));
  static const char *const nm[] = {"use", "y", "s", "bad", "odd"};
  SEXP nms = PROTECT(Rf_allocVector(STRSXP, 5));
  for (int k = 0; k < 5; ++k) SET_STRING_ELT(nms, k, Rf_mkChar(nm[k]));
  Rf_setAttrib(o, R_NamesSymbol, nms);
  UNPROTECT(7);
  return o;
}

// mcartst(): Little's test of the columns and rows of x with a value: out
// (8), mu, sigma, the columns kept (from 1) and the note; 1 (an infinite
// value), 2 (too few columns or rows) or 3 (collinear columns) instead.
SEXP C_clab_r_mcar2(SEXP x, SEXP maxit, SEXP tol) {
  clab::who = "mcartst";
  const int n = Rf_nrows(x), p = Rf_ncols(x);
  int64_t m = 0;
  int k = 0, why = 0, rc = CLAB_OK;
  SEXP r = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(r, 0, Rf_allocVector(REALSXP, 8));
  double *o = REAL(VECTOR_ELT(r, 0));
  for (int i = 0; i < 8; ++i) o[i] = 0;
  std::string note;
  {
    std::vector<double> y(static_cast<size_t>(n) * static_cast<size_t>(p), 0.0);
    std::vector<int> keep(static_cast<size_t>(p) + 1, 0);
    rc = clab_mcdrop(REAL(x), n, p, y.data(), &m, keep.data(), &k, &why);
    if (rc == CLAB_OK) {
      SET_VECTOR_ELT(r, 1, Rf_allocVector(REALSXP, k));
      SET_VECTOR_ELT(r, 2, Rf_allocMatrix(REALSXP, k, k));
      SET_VECTOR_ELT(r, 3, Rf_allocVector(INTSXP, k));
      for (int j = 0; j < k; ++j) INTEGER(VECTOR_ELT(r, 3))[j] = keep[static_cast<size_t>(j)];
      rc = clab_mcar(y.data(), static_cast<int>(m), k, Rf_asInteger(maxit), Rf_asReal(tol), o,
                     REAL(VECTOR_ELT(r, 1)), REAL(VECTOR_ELT(r, 2)));
      if (rc != CLAB_OK && o[0] == 1) why = 3;
      if (rc == CLAB_OK) note = dscd::mcnote(o);
    }
  }
  if (why > 0) {
    UNPROTECT(1);
    return Rf_ScalarInteger(why);
  }
  clab::ichk(rc);
  if (rc != CLAB_OK) {
    throw std::runtime_error("the EM estimates did not converge to a usable covariance: the columns may be "
                             "collinear, or a pattern may have too few rows.");
  }
  SET_VECTOR_ELT(r, 4, note.empty() ? Rf_ScalarString(NA_STRING) : chr({note}));
  UNPROTECT(1);
  return r;
}

}  // extern "C"
