// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the policy-evaluation R entry points (pevr.h).
#include "pevr.h"

#include <algorithm>
#include <cmath>
#include <cstring>

#include "rchk.h"
#include "spin.h"

using clab::fail;

namespace pevr {

double na(double v) { return std::isnan(v) ? NA_REAL : v; }

cols::cols(int ncol) {
  v = PROTECT(Rf_allocVector(VECSXP, ncol));
  nm.reserve(z(ncol));
}

double *cols::dbl(const char *name, int n) {
  SEXP c = Rf_allocVector(REALSXP, n);
  SET_VECTOR_ELT(v, k++, c);
  nm.emplace_back(name);
  return REAL(c);
}

SEXP cols::str(const char *name, int n) {
  SEXP c = Rf_allocVector(STRSXP, n);
  SET_VECTOR_ELT(v, k++, c);
  nm.emplace_back(name);
  return c;
}

SEXP cols::frame() {
  std::vector<const char *> p;
  for (const auto &s : nm) p.push_back(s.c_str());
  SEXP w = PROTECT(Rf_allocVector(VECSXP, k));
  for (int j = 0; j < k; ++j) SET_VECTOR_ELT(w, j, VECTOR_ELT(v, j));
  SEXP f = clab::frame(w, p.data());
  UNPROTECT(2);
  return f;
}

SEXP nlist(std::initializer_list<std::pair<const char *, SEXP>> e) {
  SEXP v = PROTECT(Rf_allocVector(VECSXP, static_cast<R_xlen_t>(e.size())));
  SEXP n = PROTECT(Rf_allocVector(STRSXP, static_cast<R_xlen_t>(e.size())));
  R_xlen_t i = 0;
  for (const auto &x : e) {
    SET_VECTOR_ELT(v, i, x.second);
    SET_STRING_ELT(n, i, Rf_mkChar(x.first));
    ++i;
  }
  Rf_setAttrib(v, R_NamesSymbol, n);
  UNPROTECT(2);
  return v;
}

void adjcols(SEXP adj, int n, std::vector<term> &out) {
  if (Rf_isNull(adj)) return;
  std::vector<std::pair<std::string, SEXP>> cs;
  if (TYPEOF(adj) == VECSXP) {
    SEXP nms = Rf_getAttrib(adj, R_NamesSymbol);
    for (R_xlen_t j = 0; j < XLENGTH(adj); ++j) {
      std::string nm = !Rf_isNull(nms) && STRING_ELT(nms, j) != NA_STRING &&
                               CHAR(STRING_ELT(nms, j))[0] != '\0'
                           ? Rf_translateCharUTF8(STRING_ELT(nms, j))
                           : "adj" + std::to_string(j + 1);
      cs.emplace_back(nm, VECTOR_ELT(adj, j));
    }
  } else if (clab::isnum(adj) && Rf_isMatrix(adj)) {
    const int nr = Rf_nrows(adj), nc = Rf_ncols(adj);
    if (nr != n) fail("`adj` must have one row per row of the panel, %d.", n);
    SEXP dn = Rf_getAttrib(adj, R_DimNamesSymbol);
    SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
    SEXP v = PROTECT(clab::dbl(adj));
    for (int j = 0; j < nc; ++j) {
      term t;
      const bool named = !Rf_isNull(cn) && STRING_ELT(cn, j) != NA_STRING &&
                         CHAR(STRING_ELT(cn, j))[0] != '\0';
      t.nm = named ? Rf_translateCharUTF8(STRING_ELT(cn, j)) : "adj" + std::to_string(j + 1);
      t.x.assign(REAL(v) + z(j) * z(n), REAL(v) + z(j + 1) * z(n));
      for (double e : t.x) {
        if (!std::isfinite(e)) fail("`adj` must hold finite numbers with no missing values.");
      }
      out.push_back(t);
    }
    UNPROTECT(1);
    return;
  } else {
    cs.emplace_back("adj", adj);
  }
  for (auto &c : cs) {
    SEXP x = c.second;
    const std::string an = "adj " + c.first;
    clab::unlab(x, an.c_str());
    if (!Rf_isVectorAtomic(x) || XLENGTH(x) != n) {
      fail("`adj` column %s must be a vector with one value per row of the panel, %d.",
           c.first.c_str(), n);
    }
    if (clab::isnum(x) && !Rf_isFactor(x)) {
      SEXP v = PROTECT(clab::num(x, an.c_str(), n));
      term t;
      t.nm = c.first;
      t.x.assign(REAL(v), REAL(v) + n);
      out.push_back(t);
      UNPROTECT(1);
      continue;
    }
    csp::lab l = csp::codes(x, an.c_str(), n);
    for (int k = 1; k < l.g; ++k) {
      term t;
      t.nm = c.first + "=" + Rf_translateCharUTF8(STRING_ELT(l.lev, k));
      t.x.resize(z(n));
      for (int i = 0; i < n; ++i) t.x[z(i)] = INTEGER(l.code)[i] == k ? 1.0 : 0.0;
      out.push_back(t);
    }
    UNPROTECT(2);
  }
}

pan panel(SEXP y, SEXP unit, SEXP time, double lo, int *np) {
  pan p;
  clab::vec(y, "y");
  SEXP yv = PROTECT(clab::num(y, "y", -1, lo));
  ++*np;
  p.n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (p.n < 4) fail("`y` must have at least four rows (two units in two periods).");
  p.y = REAL(yv);
  csp::lab l = csp::codes(unit, "unit", p.n);
  *np += 2;
  p.u.assign(INTEGER(l.code), INTEGER(l.code) + p.n);
  p.nu = l.g;
  p.ulev = l.lev;
  clab::vec(time, "time");
  SEXP tv = PROTECT(clab::num(time, "time", p.n));
  ++*np;
  const double *tt = REAL(tv);
  p.tv.assign(tt, tt + p.n);
  std::sort(p.tv.begin(), p.tv.end());
  p.tv.erase(std::unique(p.tv.begin(), p.tv.end()), p.tv.end());
  p.nt = static_cast<int>(p.tv.size());
  p.t.resize(z(p.n));
  for (int i = 0; i < p.n; ++i) {
    p.t[z(i)] = static_cast<int>(std::lower_bound(p.tv.begin(), p.tv.end(), tt[i]) - p.tv.begin());
  }
  if (p.nu < 2 || p.nt < 2) fail("the panel needs two units and two periods at least.");
  std::vector<char> seen(z(p.nu) * z(p.nt), 0);
  for (int i = 0; i < p.n; ++i) {
    char &s = seen[z(p.u[z(i)]) * z(p.nt) + z(p.t[z(i)])];
    if (s) fail("`unit` and `time` repeat a pair: one row per unit and period.");
    s = 1;
  }
  return p;
}

std::vector<int> firsts(const pan &pn, const int *d) {
  // d by unit and period, -1 where the panel has no row
  std::vector<int> dw(z(pn.nu) * z(pn.nt), -1), g(z(pn.nu), pn.nt);
  for (int i = 0; i < pn.n; ++i) dw[z(pn.u[z(i)]) * z(pn.nt) + z(pn.t[z(i)])] = d[i];
  for (int a = 0; a < pn.nu; ++a) {
    for (int s = 0; s < pn.nt; ++s) {
      const int v = dw[z(a) * z(pn.nt) + z(s)];
      if (v == 1 && g[z(a)] == pn.nt) g[z(a)] = s;
      if (v == 0 && g[z(a)] < pn.nt) {
        fail("`trt` must stay 1 once a unit is treated (unit %s).",
             Rf_translateCharUTF8(STRING_ELT(pn.ulev, a)));
      }
    }
  }
  return g;
}

std::vector<double> firstcol(SEXP first, const pan &pn) {
  if (Rf_xlength(first) != pn.n || !clab::isnum(first)) {
    fail("`first` must be numbers, one per row: the period each unit is first treated.");
  }
  SEXP f = PROTECT(clab::dbl(first));
  std::vector<double> g(z(pn.nu));
  std::vector<char> seen(z(pn.nu), 0);
  for (int i = 0; i < pn.n; ++i) {
    double v = REAL(f)[i];
    if (ISNAN(v) || v == 0 || v == R_PosInf) v = R_PosInf;
    if (v == R_NegInf) fail("`first` must be a period, or NA, 0 or Inf for never.");
    const int a = pn.u[z(i)];
    if (seen[z(a)] && !(g[z(a)] == v)) fail("`first` must be constant within each unit.");
    seen[z(a)] = 1;
    g[z(a)] = v;
  }
  UNPROTECT(1);
  return g;
}

}  // namespace pevr
