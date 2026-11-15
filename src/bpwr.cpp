// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of hwpwr() (src/pwr.h): a snapshot as text, the
// figures over a list of snapshots, and the program-output parsers on
// recorded text (for the tests). The R wrapper has checked every argument.
// Nothing calls back into the R evaluator, and no C++ object is alive when
// an error is raised.
#include <cmath>
#include <string>
#include <vector>

#include "bcol.h"
#include "pwr.h"
#include "rchk.h"

using clab::fail;

namespace {

std::string s1(SEXP x) {
  return TYPEOF(x) == STRSXP && XLENGTH(x) == 1 && STRING_ELT(x, 0) != NA_STRING ? bj::str(x, 0) : "";
}

double nn(double x) { return std::isnan(x) ? NA_REAL : x; }

// a character column of n strings, put in the frame f first
void scol(bj::cols &f, const char *name, const std::vector<std::string> &v) {
  SEXP c = Rf_allocVector(STRSXP, static_cast<R_xlen_t>(v.size()));
  f.put(name, c);
  for (std::size_t i = 0; i < v.size(); ++i) SET_STRING_ELT(c, static_cast<R_xlen_t>(i), Rf_mkCharCE(v[i].c_str(), CE_UTF8));
}

}  // namespace

extern "C" {

SEXP C_clab_r_pwrs(SEXP root, SEXP nvml, SEXP smi, SEXP ioreg, SEXP skip, SEXP pid) {
  clab::who = "hwpwr";
  std::string txt;
  bool ok = true;
  try {
    pwr::Cfg c;
    c.root = s1(root);
    c.nvml = s1(nvml);
    c.smi = s1(smi);
    c.ioreg = s1(ioreg);
    if (TYPEOF(skip) == STRSXP) {
      for (R_xlen_t i = 0; i < XLENGTH(skip); ++i) c.skip.emplace_back(bj::str(skip, i));
    }
    c.pid = TYPEOF(pid) == INTSXP && XLENGTH(pid) == 1 && INTEGER(pid)[0] != NA_INTEGER ? INTEGER(pid)[0] : 0;
    txt = pwr::put(pwr::snap(c));
  } catch (...) {
    ok = false;
  }
  if (!ok) {
    std::string().swap(txt);
    fail("the reading failed: out of memory.");
  }
  return Rf_mkString(txt.c_str());
}

SEXP C_clab_r_pwrd(SEXP snaps, SEXP rated) {
  clab::who = "hwpwr";
  std::vector<pwr::Row> r;
  double T = 0;
  int err = 0;
  try {
    std::vector<pwr::Snap> s(static_cast<std::size_t>(XLENGTH(snaps)));
    for (R_xlen_t i = 0; i < XLENGTH(snaps) && !err; ++i) {
      if (STRING_ELT(snaps, i) == NA_STRING || !pwr::get(bj::str(snaps, i), s[static_cast<std::size_t>(i)])) err = 1;
      if (!err && i > 0 && s[static_cast<std::size_t>(i)].t < s[static_cast<std::size_t>(i - 1)].t) err = 2;
    }
    if (!err && !s.empty()) {
      pwr::Rated rt;
      const double *rv = REAL(rated);
      rt.cpu = ISNAN(rv[0]) ? -1 : rv[0];
      rt.gpu = ISNAN(rv[1]) ? -1 : rv[1];
      r = pwr::rows(s, rt);
      T = s.back().t - s.front().t;
    }
  } catch (...) {
    err = 3;
  }
  if (err) std::vector<pwr::Row>().swap(r);
  if (err == 1) fail("`since` is not a reading of hwpwr().");
  if (err == 2) fail("`since` is not an earlier reading of this session: the clock went back.");
  if (err == 3) fail("the reading failed: out of memory.");
  const int n = static_cast<int>(r.size());
  std::vector<std::string> kind, ven, name, scope, how, why;
  for (const pwr::Row &x : r) {
    kind.push_back(x.kind);
    ven.push_back(x.vendor);
    name.push_back(x.name);
    scope.push_back(x.scope);
    how.push_back(x.how);
    why.push_back(x.why);
  }
  bj::cols f(11, n);
  scol(f, "kind", kind);
  scol(f, "vendor", ven);
  scol(f, "name", name);
  scol(f, "scope", scope);
  double *w = f.dbl("watts"), *j = f.dbl("joules"), *u = f.dbl("util");
  scol(f, "how", how);
  double *lo = f.dbl("flo"), *hi = f.dbl("fhi");
  scol(f, "why", why);
  for (int i = 0; i < n; ++i) {
    w[i] = nn(r[static_cast<std::size_t>(i)].watts);
    j[i] = nn(r[static_cast<std::size_t>(i)].joules);
    u[i] = nn(r[static_cast<std::size_t>(i)].util);
    lo[i] = nn(r[static_cast<std::size_t>(i)].flo);
    hi[i] = nn(r[static_cast<std::size_t>(i)].fhi);
  }
  SEXP out = PROTECT(f.done());
  SEXP dt = PROTECT(Rf_ScalarReal(T));
  Rf_setAttrib(out, Rf_install("dt"), dt);
  UNPROTECT(3);
  return out;
}

// kind "smi": text -> data frame name, watts, util, why; "ioacc": text ->
// busy shares; "iobat": text -> list(watts, why); "emi": a raw block and
// its version -> list(name, pwh).
SEXP C_clab_r_pwrp(SEXP kind, SEXP x, SEXP ver) {
  clab::who = "hwpwr";
  const std::string k = s1(kind);
  if (k == "smi" || k == "ioacc" || k == "iobat") {
    const std::string t = s1(x);
    if (k == "ioacc") {
      const std::vector<double> v = pwr::ioacc(t);
      SEXP o = PROTECT(Rf_allocVector(REALSXP, static_cast<R_xlen_t>(v.size())));
      for (std::size_t i = 0; i < v.size(); ++i) REAL(o)[i] = v[i];
      UNPROTECT(1);
      return o;
    }
    if (k == "iobat") {
      std::string why;
      const double w = pwr::iobat(t, why);
      SEXP o = PROTECT(Rf_allocVector(VECSXP, 2));
      SET_VECTOR_ELT(o, 0, Rf_ScalarReal(w < 0 ? NA_REAL : w));
      SET_VECTOR_ELT(o, 1, Rf_mkString(why.c_str()));
      static const char *const nm[] = {"watts", "why"};
      SEXP r = clab::nlist(o, nm);
      UNPROTECT(1);
      return r;
    }
    pwr::Snap s;
    pwr::smi(t, s);
    const int n = static_cast<int>(s.d.size());
    std::vector<std::string> name, why;
    for (const pwr::Dev &d : s.d) {
      name.push_back(d.name);
      why.push_back(d.why);
    }
    bj::cols f(4, n);
    scol(f, "name", name);
    double *w = f.dbl("watts"), *u = f.dbl("util");
    scol(f, "why", why);
    for (int i = 0; i < n; ++i) {
      w[i] = s.d[static_cast<std::size_t>(i)].w < 0 ? NA_REAL : s.d[static_cast<std::size_t>(i)].w;
      u[i] = s.d[static_cast<std::size_t>(i)].util < 0 ? NA_REAL : s.d[static_cast<std::size_t>(i)].util;
    }
    SEXP out = PROTECT(f.done());
    UNPROTECT(2);
    return out;
  }
  if (k == "emi" && TYPEOF(x) == RAWSXP) {
    std::vector<std::string> nm;
    std::vector<int> pwh;
    const int v = Rf_asInteger(ver);
    const bool ok = pwr::emimeta(RAW(x), static_cast<unsigned long>(XLENGTH(x)), v, nm, pwh);
    SEXP o = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP a = Rf_allocVector(STRSXP, ok ? static_cast<R_xlen_t>(nm.size()) : 0);
    SET_VECTOR_ELT(o, 0, a);
    SEXP b = Rf_allocVector(LGLSXP, ok ? static_cast<R_xlen_t>(pwh.size()) : 0);
    SET_VECTOR_ELT(o, 1, b);
    for (std::size_t i = 0; ok && i < nm.size(); ++i) {
      SET_STRING_ELT(a, static_cast<R_xlen_t>(i), Rf_mkChar(nm[i].c_str()));
      LOGICAL(b)[i] = pwh[i];
    }
    static const char *const on[] = {"name", "pwh"};
    SEXP r = clab::nlist(o, on);
    UNPROTECT(1);
    return r;
  }
  fail("unknown parser.");
}

}  // extern "C"
