// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the spatial crime methods: near repeats (src/knx.cpp)
// and the rest that follow. Each reads its points through csp::xy()
// (spin.h), checks the other arguments, calls the kernels and builds the R
// result. Nothing here calls back into the R evaluator.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

// R's distribution functions under their library names (Rmath.h declares
// them so only in newer versions)
extern "C" double Rf_ppois(double, double, int, int);
extern "C" double Rf_pnorm5(double, double, double, int, int);
extern "C" double Rf_pchisq(double, double, int, int);
extern "C" double Rf_qchisq(double, double, int, int);

using clab::fail;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

// a data frame of k double columns of length n, named nm; col[j] gives
// each column's storage
SEXP dframe(int k, int n, const char *const *nm, double **col) {
  SEXP v = PROTECT(Rf_allocVector(VECSXP, k));
  for (int j = 0; j < k; ++j) {
    SET_VECTOR_ELT(v, j, real(n));
    col[j] = REAL(VECTOR_ELT(v, j));
  }
  SEXP f = clab::frame(v, nm);
  UNPROTECT(1);
  return f;
}

// band edges: positive, finite and increasing (returned protected)
SEXP edges(SEXP e, const char *nm) {
  SEXP v = PROTECT(num(e, nm, -1, 0.0, true));
  const R_xlen_t k = XLENGTH(v);
  for (R_xlen_t q = 1; q < k; ++q) {
    if (!(REAL(v)[q] > REAL(v)[q - 1])) fail("`%s` must increase.", nm);
  }
  return v;
}

const char *const xyn[] = {"x", "y"};

// a classed list
SEXP klass(SEXP x, const char *cl) {
  Rf_classgets(x, PROTECT(Rf_mkString(cl)));
  UNPROTECT(1);
  return x;
}

// the normal-reference bandwidth of one axis, MASS::bandwidth.nrd() / 4:
// 1.06 min(sd, IQR / 1.34) n^(-1/5), IQR from type-7 quantiles
double nrd(const double *v, int n) {
  std::vector<double> u(v, v + n);
  std::sort(u.begin(), u.end());
  auto q = [&](double pr) {
    const double h = (n - 1) * pr;
    const int lo = static_cast<int>(std::floor(h));
    return lo + 1 < n ? u[lo] + (h - lo) * (u[lo + 1] - u[lo]) : u[lo];
  };
  double m = 0, ss = 0;
  for (int i = 0; i < n; ++i) m += v[i];
  m /= n;
  for (int i = 0; i < n; ++i) ss += (v[i] - m) * (v[i] - m);
  const double sd = std::sqrt(ss / (n - 1));
  return 1.06 * std::min(sd, (q(0.75) - q(0.25)) / 1.34) * std::pow(n, -0.2);
}

// the numbers of a matrix (or vector) argument: finite, at least lo,
// as doubles (returned protected)
SEXP mnum(SEXP x, const char *nm, double lo = clab::ninf) {
  clab::unlab(x, nm);
  if (!clab::isnum(x)) fail("`%s` must be numeric.", nm);
  SEXP v = PROTECT(clab::dbl(x));
  for (R_xlen_t q = 0; q < XLENGTH(v); ++q) {
    if (!std::isfinite(REAL(v)[q]) || REAL(v)[q] < lo) {
      fail("`%s` must hold finite numbers%s with no missing values.", nm, lo == 0 ? " from 0" : "");
    }
  }
  return v;
}

// an element of list x by name, or R_NilValue
SEXP elt(SEXP x, const char *nm) {
  SEXP n = Rf_getAttrib(x, R_NamesSymbol);
  for (R_xlen_t j = 0; j < Rf_xlength(n); ++j) {
    if (std::strcmp(CHAR(STRING_ELT(n, j)), nm) == 0) return VECTOR_ELT(x, j);
  }
  return R_NilValue;
}

// Weights as compressed sparse rows (rows sorted by neighbour; a pair
// listed twice or a unit its own neighbour refused): spwt()'s list, an
// spdep listw, or an n x n matrix (its non-zero entries off the diagonal).
// rp, ci, w are kept in v, protected by the caller.
struct wcsr {
  SEXP v;
  const int *rp, *ci;
  const double *w;
};

wcsr wread(SEXP w, int n) {
  std::vector<int> ei, ej;
  std::vector<double> ew;
  if (Rf_inherits(w, "listw")) {
    SEXP nb = elt(w, "neighbours"), wt = elt(w, "weights");
    if (TYPEOF(nb) != VECSXP || TYPEOF(wt) != VECSXP || XLENGTH(nb) != n || XLENGTH(wt) != n) {
      fail("`w` must have one row of neighbours and weights per unit.");
    }
    for (int i = 0; i < n; ++i) {
      SEXP a = VECTOR_ELT(nb, i), b = VECTOR_ELT(wt, i);
      if (TYPEOF(a) != INTSXP) fail("`w` must have integer neighbour lists.");
      const R_xlen_t k = XLENGTH(a);
      if (k == 1 && INTEGER(a)[0] == 0) continue;
      if (TYPEOF(b) != REALSXP || XLENGTH(b) != k) fail("`w` must have one weight per neighbour.");
      for (R_xlen_t q = 0; q < k; ++q) {
        const int j = INTEGER(a)[q];
        if (j < 1 || j > n) fail("`w` names a neighbour outside 1 to %d.", n);
        ei.push_back(i);
        ej.push_back(j - 1);
        ew.push_back(REAL(b)[q]);
      }
    }
  } else if (Rf_isMatrix(w)) {
    if (Rf_nrows(w) != n || Rf_ncols(w) != n) fail("`w` must be a %d x %d matrix.", n, n);
    SEXP m = mnum(w, "w", 0.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double v = REAL(m)[i + static_cast<R_xlen_t>(j) * n];
        if (v == 0) continue;
        if (i == j) fail("a unit cannot be its own neighbour in `w`.");
        ei.push_back(i);
        ej.push_back(j);
        ew.push_back(v);
      }
    }
    UNPROTECT(1);
  } else {
    SEXP nn = TYPEOF(w) == VECSXP ? elt(w, "n") : R_NilValue;
    if (Rf_isNull(nn)) fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
    if (clab::val1(nn) != n) fail("`w` is for %g units, not %d.", clab::val1(nn), n);
    SEXP i = PROTECT(num(elt(w, "i"), "w$i", -1, 1.0));
    const R_xlen_t m = XLENGTH(i);
    SEXP j = PROTECT(num(elt(w, "j"), "w$j", m, 1.0));
    SEXP v = PROTECT(num(elt(w, "w"), "w$w", m, 0.0));
    for (R_xlen_t q = 0; q < m; ++q) {
      const double a = REAL(i)[q], b = REAL(j)[q];
      if (a > n || b > n || a != std::floor(a) || b != std::floor(b)) {
        fail("`w$i` and `w$j` must be unit numbers from 1 to %d.", n);
      }
      ei.push_back(static_cast<int>(a) - 1);
      ej.push_back(static_cast<int>(b) - 1);
      ew.push_back(REAL(v)[q]);
    }
    UNPROTECT(3);
  }
  const std::size_t m = ei.size();
  for (std::size_t q = 0; q < m; ++q) {
    if (!(ew[q] >= 0) || !std::isfinite(ew[q])) fail("`w` must hold finite weights from 0.");
  }
  std::vector<int> o(m);
  for (std::size_t q = 0; q < m; ++q) o[q] = static_cast<int>(q);
  std::sort(o.begin(), o.end(),
            [&](int a, int b) { return ei[a] < ei[b] || (ei[a] == ei[b] && ej[a] < ej[b]); });
  wcsr c{PROTECT(Rf_allocVector(VECSXP, 3)), nullptr, nullptr, nullptr};
  SET_VECTOR_ELT(c.v, 0, Rf_allocVector(INTSXP, n + 1));
  SET_VECTOR_ELT(c.v, 1, Rf_allocVector(INTSXP, static_cast<R_xlen_t>(m)));
  SET_VECTOR_ELT(c.v, 2, real(static_cast<R_xlen_t>(m)));
  int *rp = INTEGER(VECTOR_ELT(c.v, 0)), *ci = INTEGER(VECTOR_ELT(c.v, 1));
  double *ww = REAL(VECTOR_ELT(c.v, 2));
  std::fill(rp, rp + n + 1, 0);
  for (std::size_t q = 0; q < m; ++q) {
    const int a = o[q];
    if (ei[a] == ej[a]) fail("a unit cannot be its own neighbour in `w`.");
    if (q > 0 && ei[a] == ei[o[q - 1]] && ej[a] == ej[o[q - 1]])
      fail("`w` lists a pair of units twice.");
    ci[q] = ej[a];
    ww[q] = ew[a];
    ++rp[ei[a] + 1];
  }
  for (int i = 0; i < n; ++i) rp[i + 1] += rp[i];
  c.rp = rp;
  c.ci = ci;
  c.w = ww;
  return c;
}

// the names of the columns of matrix x, or prefix and number
SEXP cnames(SEXP x, int p, const char *pre) {
  SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
  SEXP cn = Rf_isNull(dn) ? R_NilValue : VECTOR_ELT(dn, 1);
  SEXP out = PROTECT(Rf_allocVector(STRSXP, p));
  char buf[32];
  for (int j = 0; j < p; ++j) {
    if (!Rf_isNull(cn)) {
      SET_STRING_ELT(out, j, STRING_ELT(cn, j));
    } else {
      std::snprintf(buf, sizeof buf, "%s%d", pre, j + 1);
      SET_STRING_ELT(out, j, Rf_mkChar(buf));
    }
  }
  UNPROTECT(1);
  return out;
}

}  // namespace

extern "C" {

// ---------------------------------------------------------- near repeats

SEXP C_clab_r_knxnr(SEXP x, SEXP y, SEXP t, SEXP ds, SEXP dt, SEXP nsim, SEXP seed, SEXP lonlat) {
  clab::who = "knxnr";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 incidents are needed.");
  SEXP tv = PROTECT(num(t, "t", n));
  SEXP sv = edges(ds, "ds");
  SEXP bv = edges(dt, "dt");
  const int ks = static_cast<int>(XLENGTH(sv)), kt = static_cast<int>(XLENGTH(bv));
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  const int m = ks * kt;
  SEXP o = PROTECT(real(3 * static_cast<R_xlen_t>(m) + 1));
  kc(clab_knxnr(p.x, p.y, REAL(tv), n, REAL(sv), ks, REAL(bv), kt, ns, sd, REAL(o)),
     "the kernel refused the input.");
  static const char *const tn[] = {"dlo", "dhi", "tlo", "thi", "obs", "exp", "ratio", "p", "ppois"};
  double *c[9];
  SEXP tab = PROTECT(dframe(9, m, tn, c));
  const double *ob = REAL(o), *ex = REAL(o) + m, *pp = REAL(o) + 2 * m;
  for (int b = 0; b < kt; ++b) {
    for (int s = 0; s < ks; ++s) {
      const int q = s + b * ks;
      c[0][q] = s ? REAL(sv)[s - 1] : 0.0;
      c[1][q] = REAL(sv)[s];
      c[2][q] = b ? REAL(bv)[b - 1] : 0.0;
      c[3][q] = REAL(bv)[b];
      c[4][q] = ob[q];
      c[5][q] = ex[q];
      c[6][q] = ex[q] > 0 ? ob[q] / ex[q] : std::numeric_limits<double>::quiet_NaN();
      c[7][q] = pp[q];
      // P(X >= obs) for X Poisson with the expected count as its mean
      c[8][q] = ob[q] > 0 ? Rf_ppois(ob[q] - 1, ex[q], 0, 0) : 1.0;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, tab);
  SET_VECTOR_ELT(out, 1, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 2, Rf_ScalarReal(REAL(o)[3 * m]));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ns));
  static const char *const on[] = {"tab", "n", "pairs", "nsim"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "knxnr");
}

// --------------------------------------------------------------- hot spots

SEXP C_clab_r_kdsurf(SEXP x, SEXP y, SEXP w, SEXP h, SEXP win, SEXP nx, SEXP ny, SEXP edge, SEXP at,
                     SEXP cover) {
  clab::who = "kdsurf";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 2) fail("at least 2 points are needed.");
  SEXP wv = PROTECT(Rf_isNull(w) ? R_NilValue : num(w, "w", n, 0.0));
  double hh[2];
  if (Rf_isNull(h)) {
    hh[0] = nrd(p.x, n);
    hh[1] = nrd(p.y, n);
    if (!(hh[0] > 0 && hh[1] > 0)) {
      fail("the points do not spread along both axes: give `h`.");
    }
  } else {
    SEXP hv = PROTECT(num(h, "h", -1, 0.0, true));
    if (XLENGTH(hv) != 1 && XLENGTH(hv) != 2) fail("`h` must be one or two positive numbers.");
    hh[0] = REAL(hv)[0];
    hh[1] = REAL(hv)[XLENGTH(hv) - 1];
    UNPROTECT(1);
  }
  static const char *const eo[] = {"none", "uniform", "diggle"};
  const int ed = csp::pick(edge, "edge", eo, 3);
  double wn[4];
  const bool hw = !Rf_isNull(win);
  if (hw) csp::win(win, wn);
  if (ed > 0 && !hw) fail("an edge correction needs `win`.");
  const double *pw = Rf_isNull(wv) ? nullptr : REAL(wv);
  SEXP out;
  if (!Rf_isNull(at)) {
    const bool own = TYPEOF(at) == STRSXP;
    if (own && !(XLENGTH(at) == 1 && std::strcmp(CHAR(STRING_ELT(at, 0)), "points") == 0)) {
      fail("`at` must be NULL, \"points\" or coordinates.");
    }
    csp::pts q = own ? p : csp::xy(at, R_NilValue, R_NilValue, csp::planar);
    SEXP z = PROTECT(real(q.n));
    kc(clab_kdsurf(p.x, p.y, pw, n, hh[0], hh[1], q.x, q.n, q.y, q.n, 0, hw ? wn : nullptr, ed,
                   own ? 1 : 0, REAL(z)),
       "a point lies where its window holds none of its kernel.");
    out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, clab::dcopy(VECTOR_ELT(q.v, 0)));
    SET_VECTOR_ELT(out, 1, clab::dcopy(VECTOR_ELT(q.v, 1)));
    SET_VECTOR_ELT(out, 2, z);
    SET_VECTOR_ELT(out, 3, clab::nvec(2, xyn, hh));
    SET_VECTOR_ELT(out, 4, R_NilValue);
    UNPROTECT(own ? 2 : 3);
  } else {
    const int ax = clab::whole(nx, "nx", 2), ay = clab::whole(ny, "ny", 2);
    if (static_cast<double>(ax) * ay > 1e8) fail("`nx` times `ny` must be at most 1e8.");
    SEXP cv = PROTECT(num(cover, "cover", -1, 0.0, true));
    for (R_xlen_t q = 0; q < XLENGTH(cv); ++q) {
      if (!(REAL(cv)[q] < 1)) fail("`cover` must hold shares between 0 and 1.");
    }
    double b[4];
    if (hw) {
      std::memcpy(b, wn, sizeof b);
    } else {
      b[0] = *std::min_element(p.x, p.x + n) - 3 * hh[0];
      b[1] = *std::max_element(p.x, p.x + n) + 3 * hh[0];
      b[2] = *std::min_element(p.y, p.y + n) - 3 * hh[1];
      b[3] = *std::max_element(p.y, p.y + n) + 3 * hh[1];
    }
    SEXP gx = PROTECT(real(ax)), gy = PROTECT(real(ay));
    for (int a = 0; a < ax; ++a) REAL(gx)[a] = b[0] + a * (b[1] - b[0]) / (ax - 1);
    for (int a = 0; a < ay; ++a) REAL(gy)[a] = b[2] + a * (b[3] - b[2]) / (ay - 1);
    SEXP z = PROTECT(Rf_allocMatrix(REALSXP, ax, ay));
    kc(clab_kdsurf(p.x, p.y, pw, n, hh[0], hh[1], REAL(gx), ax, REAL(gy), ay, 1, hw ? wn : nullptr,
                   ed, 0, REAL(z)),
       "a point lies where its window holds none of its kernel.");
    // the smallest sets of cells holding each share of the surface's mass
    const R_xlen_t m = static_cast<R_xlen_t>(ax) * ay;
    std::vector<double> zs(REAL(z), REAL(z) + m);
    std::sort(zs.begin(), zs.end(), [](double a, double c) { return a > c; });
    double tot = 0;
    for (R_xlen_t q = 0; q < m; ++q) tot += zs[q];
    const double ca = (b[1] - b[0]) / (ax - 1) * (b[3] - b[2]) / (ay - 1);
    const int nc = static_cast<int>(XLENGTH(cv));
    static const char *const rn[] = {"cover", "level", "area", "share"};
    double *rc[4];
    SEXP reg = PROTECT(dframe(4, nc, rn, rc));
    for (int k = 0; k < nc; ++k) {
      double run = 0;
      R_xlen_t q = 0;
      while (q < m && run < REAL(cv)[k] * tot) run += zs[q++];
      rc[0][k] = REAL(cv)[k];
      rc[1][k] = zs[q > 0 ? q - 1 : 0];
      rc[2][k] = q * ca;
      rc[3][k] = static_cast<double>(q) / m;
    }
    out = PROTECT(Rf_allocVector(VECSXP, 5));
    SET_VECTOR_ELT(out, 0, gx);
    SET_VECTOR_ELT(out, 1, gy);
    SET_VECTOR_ELT(out, 2, z);
    SET_VECTOR_ELT(out, 3, clab::nvec(2, xyn, hh));
    SET_VECTOR_ELT(out, 4, reg);
    UNPROTECT(6);
  }
  PROTECT(out);
  static const char *const on[] = {"x", "y", "z", "h", "regions"};
  clab::nlist(out, on);
  UNPROTECT(3);
  return klass(out, "kdsurf");
}

SEXP C_clab_r_dbscn(SEXP x, SEXP y, SEXP eps, SEXP minpts, SEXP border, SEXP lonlat) {
  clab::who = "dbscn";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  const double e = clab::one(eps, "eps", 0.0, std::numeric_limits<double>::infinity());
  const int mp = clab::whole(minpts, "minpts", 1);
  const bool bd = csp::flag(border, "border");
  SEXP cl = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP co = PROTECT(Rf_allocVector(LGLSXP, n));
  int k = 0;
  kc(clab_dbscn(p.x, p.y, n, e, mp, bd ? 1 : 0, INTEGER(cl), LOGICAL(co), &k),
     "the kernel refused the input.");
  SEXP sz = PROTECT(Rf_allocVector(INTSXP, k));
  std::fill(INTEGER(sz), INTEGER(sz) + k, 0);
  for (int i = 0; i < n; ++i) {
    if (INTEGER(cl)[i] >= 0) {
      ++INTEGER(sz)[INTEGER(cl)[i]];
      ++INTEGER(cl)[i];
    } else {
      INTEGER(cl)[i] = NA_INTEGER;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, co);
  SET_VECTOR_ELT(out, 2, sz);
  static const char *const on[] = {"cl", "core", "sizes"};
  clab::nlist(out, on);
  UNPROTECT(5);
  return klass(out, "dbscn");
}

SEXP C_clab_r_hsper(SEXP x, SEXP y, SEXP prd, SEXP cell, SEXP top, SEXP org, SEXP lonlat) {
  clab::who = "hsper";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  csp::lab g = csp::codes(prd, "prd", n);
  const double cs = clab::one(cell, "cell", 0.0, std::numeric_limits<double>::infinity());
  const double tp = clab::val1(top);
  if (!(tp > 0 && tp <= 1)) fail("`top` must be one share above 0, at most 1.");
  double o[2];
  if (Rf_isNull(org)) {
    o[0] = *std::min_element(p.x, p.x + n);
    o[1] = *std::min_element(p.y, p.y + n);
  } else {
    SEXP ov = PROTECT(num(org, "org", 2));
    o[0] = REAL(ov)[0];
    o[1] = REAL(ov)[1];
    UNPROTECT(1);
  }
  const int np = g.g;
  SEXP cid = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP ix = PROTECT(real(n)), iy = PROTECT(real(n));
  SEXP hot = PROTECT(Rf_allocVector(INTSXP, static_cast<R_xlen_t>(n) * np));
  SEXP per = PROTECT(real(9 * static_cast<R_xlen_t>(np)));
  SEXP cv = PROTECT(real(3 * static_cast<R_xlen_t>(n)));
  int nc = 0;
  kc(clab_hsper(p.x, p.y, INTEGER(g.code), n, np, cs, o[0], o[1], tp, INTEGER(cid), REAL(ix),
                REAL(iy), &nc, INTEGER(hot), REAL(per), REAL(cv)),
     "the cells are too many for the grid: choose a larger `cell` or an `org` nearer the points.");
  static const char *const pn[] = {"period", "n",      "active",  "nhalf",  "shalf",
                                   "hot",    "tiedin", "tiedout", "hotshr", "jacc"};
  SEXP pc = PROTECT(Rf_allocVector(VECSXP, 10));
  SET_VECTOR_ELT(pc, 0, g.lev);
  for (int j = 0; j < 9; ++j) {
    SEXP c = real(np);
    SET_VECTOR_ELT(pc, j + 1, c);
    for (int q = 0; q < np; ++q) REAL(c)[q] = REAL(per)[9 * static_cast<R_xlen_t>(q) + j];
  }
  SEXP pf = PROTECT(clab::frame(pc, pn));
  static const char *const cn[] = {"x", "y", "n", "nhot", "run"};
  double *cc[5];
  SEXP cf = PROTECT(dframe(5, nc, cn, cc));
  for (int j = 0; j < nc; ++j) {
    cc[0][j] = o[0] + (REAL(ix)[j] + 0.5) * cs;
    cc[1][j] = o[1] + (REAL(iy)[j] + 0.5) * cs;
    for (int q = 0; q < 3; ++q) cc[2 + q][j] = REAL(cv)[3 * static_cast<R_xlen_t>(j) + q];
  }
  csp::back(p, cc[0], cc[1], nc);
  SEXP hm = PROTECT(Rf_allocMatrix(LGLSXP, nc, np));
  for (int q = 0; q < np; ++q) {
    for (int j = 0; j < nc; ++j) {
      LOGICAL(hm)
      [j + static_cast<R_xlen_t>(q) * nc] = INTEGER(hot)[j + static_cast<R_xlen_t>(q) * nc];
    }
  }
  for (int i = 0; i < n; ++i) ++INTEGER(cid)[i];
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, pf);
  SET_VECTOR_ELT(out, 1, cf);
  SET_VECTOR_ELT(out, 2, hm);
  SET_VECTOR_ELT(out, 3, cid);
  SET_VECTOR_ELT(out, 4, Rf_ScalarInteger(std::max(1, static_cast<int>(std::ceil(tp * nc)))));
  static const char *const on[] = {"periods", "cells", "hot", "cell", "k"};
  clab::nlist(out, on);
  UNPROTECT(14);
  return klass(out, "hsper");
}

SEXP C_clab_r_ebshr(SEXP y, SEXP v, SEXP w, SEXP b) {
  clab::who = "ebshr";
  SEXP yv = PROTECT(num(y, "y"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (n < 2) fail("at least 2 places are needed.");
  SEXP wv = PROTECT(Rf_isNull(w) ? R_NilValue : num(w, "w", n, 0.0));
  const double *pw = Rf_isNull(wv) ? nullptr : REAL(wv);
  double nv;
  if (Rf_isNull(v)) {
    double sw = 0, sy = 0;
    for (int i = 0; i < n; ++i) {
      sw += pw ? pw[i] : 1.0;
      sy += (pw ? pw[i] : 1.0) * REAL(yv)[i];
    }
    if (!(sw > 0)) fail("`w` must have a positive sum.");
    nv = sy / sw;
    if (!(nv >= 0)) fail("the mean of `y` is negative: give the noise variance `v`.");
  } else {
    SEXP vv = PROTECT(num(v, "v", 1, 0.0));
    nv = REAL(vv)[0];
    UNPROTECT(1);
  }
  const double bb = Rf_isNull(b) ? std::numeric_limits<double>::quiet_NaN() : clab::p01(b, "b");
  double fit[6];
  SEXP sh = PROTECT(real(n)), fa = PROTECT(real(n));
  kc(clab_ebshr(REAL(yv), pw, n, nv, bb, fit, REAL(sh), REAL(fa)), "`w` must have a positive sum.");
  static const char *const fn[] = {"mean", "total", "noise", "signal", "bstar", "b"};
  static const char *const en[] = {"y", "shrunk", "fall"};
  SEXP ec = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(ec, 0, clab::dcopy(yv));
  SET_VECTOR_ELT(ec, 1, sh);
  SET_VECTOR_ELT(ec, 2, fa);
  SEXP est = PROTECT(clab::frame(ec, en));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, clab::nvec(6, fn, fit));
  SET_VECTOR_ELT(out, 1, est);
  static const char *const on[] = {"fit", "est"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "ebshr");
}

// ------------------------------------------------------------------ scans

SEXP C_clab_r_kscan(SEXP x, SEXP y, SEXP t, SEXP r, SEXP tw, SEXP minc, SEXP top, SEXP nsim,
                    SEXP seed, SEXP lonlat) {
  clab::who = "kscan";
  csp::pts p = csp::xy(x, y, lonlat);
  const int n = p.n;
  if (n < 3) fail("at least 3 incidents are needed.");
  SEXP tv = PROTECT(num(t, "t", n));
  SEXP rv = edges(r, "r");
  SEXP wv = PROTECT(num(tw, "tw", -1, 0.0, true));
  const int mc = clab::whole(minc, "minc", 1);
  const int tp = clab::whole(top, "top", 1);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  const int nr = static_cast<int>(XLENGTH(rv)), nw = static_cast<int>(XLENGTH(wv));
  const int tq = tp < n ? tp : n;
  SEXP o = PROTECT(real(10 * static_cast<R_xlen_t>(tq)));
  int nc = 0;
  kc(clab_kscan(p.x, p.y, REAL(tv), n, REAL(rv), nr, REAL(wv), nw, mc, ns, sd, tq, REAL(o), &nc),
     "the kernel refused the input.");
  static const char *const cn[] = {"x",   "y",     "r",    "start", "end", "n",
                                   "exp", "ncirc", "nwin", "rr",    "llr", "p"};
  double *c[12];
  SEXP cl = PROTECT(dframe(12, nc, cn, c));
  SEXP mb = PROTECT(Rf_allocVector(INTSXP, n));
  for (int i = 0; i < n; ++i) INTEGER(mb)[i] = NA_INTEGER;
  for (int k = nc - 1; k >= 0; --k) {
    const double *v = REAL(o) + 10 * k;
    const int ci = static_cast<int>(v[0]);
    c[0][k] = p.x[ci];
    c[1][k] = p.y[ci];
    for (int j = 0; j < 7; ++j) c[2 + j][k] = v[1 + j];
    // the relative risk inside against outside the cylinder
    c[9][k] = (v[4] / v[5]) / ((n - v[4]) / (n - v[5]));
    c[10][k] = v[8];
    c[11][k] = v[9];
    for (int i = 0; i < n; ++i) {
      const double dx = p.x[i] - p.x[ci], dy = p.y[i] - p.y[ci];
      if (dx * dx + dy * dy <= v[1] * v[1] && REAL(tv)[i] >= v[2] && REAL(tv)[i] < v[3]) {
        INTEGER(mb)[i] = k + 1;
      }
    }
  }
  csp::back(p, c[0], c[1], nc);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, cl);
  SET_VECTOR_ELT(out, 1, mb);
  SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 3, Rf_ScalarInteger(ns));
  static const char *const on[] = {"cl", "memb", "n", "nsim"};
  clab::nlist(out, on);
  UNPROTECT(8);
  return klass(out, "kscan");
}

// ---------------------------------------------- spatial econometrics

SEXP C_clab_r_gmoran(SEXP x, SEXP w, SEXP nsim, SEXP seed) {
  clab::who = "gmoran";
  const bool mat = Rf_isMatrix(x);
  SEXP xv = mat ? mnum(x, "x") : PROTECT(num(x, "x"));
  const int n = mat ? Rf_nrows(x) : clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  const int k = mat ? Rf_ncols(x) : 1;
  if (n < 4) fail("at least 4 units are needed.");
  wcsr c = wread(w, n);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  SEXP o = PROTECT(real(10 * static_cast<R_xlen_t>(k)));
  kc(clab_gmoran(REAL(xv), n, k, c.rp, c.ci, c.w, ns, sd, REAL(o)),
     "`w` has no positive weight: Moran's I is undefined.");
  static const char *const nm[] = {"var", "i",  "e",  "vn",    "vr", "zn",
                                   "zr",  "pn", "pr", "pperm", "s0"};
  SEXP cols = PROTECT(Rf_allocVector(VECSXP, 11));
  SET_VECTOR_ELT(cols, 0, mat ? cnames(x, k, "x") : Rf_mkString("x"));
  for (int j = 0; j < 10; ++j) {
    SEXP v = real(k);
    SET_VECTOR_ELT(cols, j + 1, v);
    for (int q = 0; q < k; ++q) REAL(v)[q] = REAL(o)[10 * static_cast<R_xlen_t>(q) + j];
  }
  SEXP f = clab::frame(cols, nm);
  UNPROTECT(4);
  return f;
}

SEXP C_clab_r_bvmor(SEXP x, SEXP y, SEXP w, SEXP nsim, SEXP seed) {
  clab::who = "bvmor";
  SEXP xv = PROTECT(num(x, "x"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  if (n < 3) fail("at least 3 units are needed.");
  SEXP yv = PROTECT(num(y, "y", n));
  wcsr c = wread(w, n);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int sd = ns > 0 ? csp::seed(seed) : 0u;
  double o[6];
  kc(clab_bvmor(REAL(xv), REAL(yv), n, c.rp, c.ci, c.w, ns, sd, o),
     "`x` or `y` is constant, or `w` has no positive weight: the statistic is undefined.");
  static const char *const nm[] = {"i", "mean", "sd", "z", "p", "ptwo"};
  UNPROTECT(3);
  return clab::nvec(6, nm, o);
}

SEXP C_clab_r_sarml(SEXP y, SEXP x, SEXP w, SEXP model, SEXP durbin, SEXP method, SEXP icpt) {
  clab::who = "sarml";
  SEXP yv = PROTECT(num(y, "y"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(yv)), "y");
  if (n < 4) fail("at least 4 units are needed.");
  if (n > 5000) fail("at most 5000 units: the likelihood takes a dense determinant.");
  if (!Rf_isMatrix(x) || Rf_nrows(x) != n) fail("`x` must have one row per value of `y`.");
  SEXP xv = mnum(x, "x");
  const int px = Rf_ncols(x);
  SEXP xn = PROTECT(cnames(x, px, "x"));
  wcsr c = wread(w, n);
  static const char *const mo[] = {"lag", "error", "sac"};
  static const char *const me[] = {"ml", "2sls"};
  const int md = csp::pick(model, "model", mo, 3);
  const int mt = csp::pick(method, "method", me, 2);
  const bool db = csp::flag(durbin, "durbin");
  const bool ic = csp::flag(icpt, "int");
  if (mt == 1 && md != 0) fail("two-stage least squares fits the lag model only.");
  // the lags W x of the columns that vary
  std::vector<int> var;
  for (int j = 0; j < px; ++j) {
    const double *v = REAL(xv) + static_cast<R_xlen_t>(j) * n;
    bool same = true;
    for (int i = 1; i < n && same; ++i) same = v[i] == v[0];
    if (!same) var.push_back(j);
  }
  const int nv = static_cast<int>(var.size());
  auto lag = [&](const double *v, double *o) {
    for (int i = 0; i < n; ++i) {
      double s = 0;
      for (int q = c.rp[i]; q < c.rp[i + 1]; ++q) s += c.w[q] * v[c.ci[q]];
      o[i] = s;
    }
  };
  const int p = (ic ? 1 : 0) + px + (db ? nv : 0);
  if (p >= n) fail("`x` has too many columns for the number of units.");
  SEXP X = PROTECT(real(static_cast<R_xlen_t>(n) * p));
  SEXP tn = PROTECT(Rf_allocVector(STRSXP, p + 2));
  int col = 0;
  if (ic) {
    std::fill(REAL(X), REAL(X) + n, 1.0);
    SET_STRING_ELT(tn, col++, Rf_mkChar("(Intercept)"));
  }
  for (int j = 0; j < px; ++j) {
    std::memcpy(REAL(X) + static_cast<R_xlen_t>(col) * n, REAL(xv) + static_cast<R_xlen_t>(j) * n,
                n * sizeof(double));
    SET_STRING_ELT(tn, col++, STRING_ELT(xn, j));
  }
  if (db) {
    for (int j : var) {
      lag(REAL(xv) + static_cast<R_xlen_t>(j) * n, REAL(X) + static_cast<R_xlen_t>(col) * n);
      const char *nm = CHAR(STRING_ELT(xn, j));
      std::vector<char> b(std::strlen(nm) + 5);
      std::snprintf(b.data(), b.size(), "lag.%s", nm);
      SET_STRING_ELT(tn, col++, Rf_mkChar(b.data()));
    }
  }
  SET_STRING_ELT(tn, p, Rf_mkChar("rho"));
  SET_STRING_ELT(tn, p + 1, Rf_mkChar("lambda"));
  // instruments of two-stage least squares: W x and W^2 x of the varying
  // columns (without the lags already in X)
  const int ni = mt == 1 ? (db ? nv : 2 * nv) : 0;
  SEXP H = PROTECT(real(static_cast<R_xlen_t>(n) * (ni > 0 ? ni : 1)));
  if (mt == 1) {
    if (nv == 0) fail("two-stage least squares needs a column of `x` that varies.");
    std::vector<double> t(n);
    int h = 0;
    for (int j : var) {
      lag(REAL(xv) + static_cast<R_xlen_t>(j) * n, t.data());
      if (!db) std::memcpy(REAL(H) + static_cast<R_xlen_t>(h++) * n, t.data(), n * sizeof(double));
      lag(t.data(), REAL(H) + static_cast<R_xlen_t>(h++) * n);
    }
  }
  const int k = p + 2;
  SEXP cf = PROTECT(real(4 * static_cast<R_xlen_t>(p) + 4));
  SEXP vc = PROTECT(Rf_allocMatrix(REALSXP, k, k));
  double st[5];
  kc(clab_sarml(REAL(yv), REAL(X), n, p, c.rp, c.ci, c.w, md, mt, mt == 1 ? REAL(H) : nullptr, ni,
                REAL(cf), REAL(vc), st),
     "the model cannot be fitted: `x` is collinear (after filtering), or `w` has no positive "
     "weight.");
  // the coefficient table: beta, then rho and/or lambda
  const int nk = p + (md == 2 ? 2 : 1);
  static const char *const cn[] = {"term", "est", "se", "z", "p", "ols", "olsse"};
  SEXP cols = PROTECT(Rf_allocVector(VECSXP, 7));
  SEXP term = Rf_allocVector(STRSXP, nk);
  SET_VECTOR_ELT(cols, 0, term);
  double *cc[6];
  for (int j = 0; j < 6; ++j) {
    SET_VECTOR_ELT(cols, j + 1, real(nk));
    cc[j] = REAL(VECTOR_ELT(cols, j + 1));
  }
  const double *cv = REAL(cf);
  int r = 0;
  for (int j = 0; j < k; ++j) {
    if ((j == p && md == 1) || (j == p + 1 && md == 0)) continue;
    SET_STRING_ELT(term, r, STRING_ELT(tn, j));
    cc[0][r] = cv[j];
    cc[1][r] = cv[k + j];
    cc[2][r] = cv[j] / cv[k + j];
    cc[3][r] = 2 * Rf_pnorm5(-std::fabs(cc[2][r]), 0, 1, 1, 0);
    cc[4][r] = j < p ? cv[2 * p + 4 + j] : std::numeric_limits<double>::quiet_NaN();
    cc[5][r] = j < p ? cv[3 * p + 4 + j] : std::numeric_limits<double>::quiet_NaN();
    ++r;
  }
  SEXP tab = PROTECT(clab::frame(cols, cn));
  Rf_setAttrib(vc, R_DimNamesSymbol, PROTECT(Rf_allocVector(VECSXP, 2)));
  SET_VECTOR_ELT(Rf_getAttrib(vc, R_DimNamesSymbol), 0, tn);
  SET_VECTOR_ELT(Rf_getAttrib(vc, R_DimNamesSymbol), 1, tn);
  static const char *const fn[] = {"loglik", "aic", "s2", "llols", "n"};
  const double fv[5] = {st[0], st[2], st[1], st[3], static_cast<double>(n)};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 8));
  SET_VECTOR_ELT(out, 0, tab);
  SET_VECTOR_ELT(out, 1, clab::nvec(5, fn, fv));
  SET_VECTOR_ELT(out, 2, vc);
  SET_VECTOR_ELT(out, 3, Rf_mkString(mo[md]));
  SET_VECTOR_ELT(out, 4, Rf_mkString(me[mt]));
  SET_VECTOR_ELT(out, 5, Rf_ScalarLogical(db));
  SET_VECTOR_ELT(out, 6, Rf_ScalarInteger(db ? nv : 0));
  SET_VECTOR_ELT(out, 7, w);
  static const char *const on[] = {"coef", "fit", "vcov", "model", "method", "durbin", "nlag", "w"};
  clab::nlist(out, on);
  UNPROTECT(13);
  return klass(out, "sarml");
}

SEXP C_clab_r_spill(SEXP w, SEXP rho, SEXP b, SEXP th, SEXP vc, SEXP nm) {
  clab::who = "spill";
  SEXP bv = PROTECT(num(b, "b"));
  const int k = static_cast<int>(XLENGTH(bv));
  SEXP tv = PROTECT(Rf_isNull(th) ? R_NilValue : num(th, "th", k));
  const double r = clab::one(rho, "rho", -HUGE_VAL, HUGE_VAL);
  // the units: from the weights
  int n = 0;
  if (Rf_inherits(w, "listw")) {
    n = static_cast<int>(XLENGTH(elt(w, "neighbours")));
  } else if (Rf_isMatrix(w)) {
    n = Rf_nrows(w);
  } else if (TYPEOF(w) == VECSXP && !Rf_isNull(elt(w, "n"))) {
    n = static_cast<int>(clab::val1(elt(w, "n")));
  } else {
    fail("`w` must be the list spwt() returns, an spdep listw or a matrix.");
  }
  if (n < 2) fail("`w` must be for at least 2 units.");
  if (n > 5000) fail("at most 5000 units: the effects take a dense inverse.");
  wcsr c = wread(w, n);
  // vc comes from spill(), which takes it from the fit, m x m
  SEXP vv = Rf_isNull(vc) ? PROTECT(R_NilValue) : mnum(vc, "vcov");
  SEXP o = PROTECT(real(6 * static_cast<R_xlen_t>(k)));
  kc(clab_spill(n, c.rp, c.ci, c.w, r, REAL(bv), Rf_isNull(tv) ? nullptr : REAL(tv), k,
                Rf_isNull(vv) ? nullptr : REAL(vv), REAL(o)),
     "I - rho W is singular: no spillover multiplier.");
  static const char *const cn[] = {"term", "direct", "indirect", "total", "sed", "sei", "set"};
  SEXP cols = PROTECT(Rf_allocVector(VECSXP, 7));
  SET_VECTOR_ELT(cols, 0, Rf_isNull(nm) ? cnames(R_NilValue, k, "x") : nm);
  for (int j = 0; j < 6; ++j) {
    SEXP v = real(k);
    SET_VECTOR_ELT(cols, j + 1, v);
    for (int q = 0; q < k; ++q) REAL(v)[q] = REAL(o)[6 * static_cast<R_xlen_t>(q) + j];
  }
  SEXP f = clab::frame(cols, cn);
  UNPROTECT(6);
  return f;
}

// ---------------------------------------------------- point patterns

SEXP C_clab_r_crossk(SEXP x, SEXP y, SEXP marks, SEXP i, SEXP j, SEXP win, SEXP r, SEXP sigma,
                     SEXP nsim, SEXP seed) {
  clab::who = "crossk";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (n < 2) fail("at least 2 points are needed.");
  csp::lab g = csp::codes(marks, "marks", n);
  const int ti = clab::ref(i, g.lev, g.g, "i", "type");
  const int tj = clab::ref(j, g.lev, g.g, "j", "type");
  double wn[4];
  csp::win(win, wn);
  for (int a = 0; a < n; ++a) {
    if (!csp::inside(wn, p.x[a], p.y[a])) fail("every point must lie inside `win`.");
  }
  SEXP rv = PROTECT(num(r, "r", -1, 0.0));
  for (R_xlen_t q = 1; q < XLENGTH(rv); ++q) {
    if (!(REAL(rv)[q] > REAL(rv)[q - 1])) fail("`r` must increase.");
  }
  const int nr = static_cast<int>(XLENGTH(rv));
  const bool sg = !Rf_isNull(sigma);
  const double sd = sg ? clab::one(sigma, "sigma", 0.0, HUGE_VAL) : 0.0;
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int se = sg && ns > 0 ? csp::seed(seed) : 0u;
  static const char *const kn[] = {"r", "k", "l", "theo"};
  double *c[4];
  SEXP kf = PROTECT(dframe(4, nr, kn, c));
  double seg[2];
  kc(clab_crossk(p.x, p.y, INTEGER(g.code), n, ti - 1, tj - 1, wn, REAL(rv), nr, sd, ns, se, c[1],
                 seg),
     "type `i` needs a point and type `j` another, or the kernel weights of a point vanish: "
     "widen `sigma`.");
  for (int q = 0; q < nr; ++q) {
    const double rr = REAL(rv)[q];
    c[0][q] = rr;
    c[2][q] = std::sqrt(c[1][q] / 3.14159265358979323846);
    c[3][q] = 3.14159265358979323846 * rr * rr;
  }
  SEXP sv = R_NilValue;
  if (sg) {
    static const char *const sn[] = {"t", "p", "nsim"};
    const double v[3] = {seg[0], seg[1], static_cast<double>(ns)};
    sv = clab::nvec(3, sn, v);
  }
  PROTECT(sv);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, kf);
  SET_VECTOR_ELT(out, 1, sv);
  SET_VECTOR_ELT(out, 2, Rf_ScalarString(STRING_ELT(g.lev, ti - 1)));
  SET_VECTOR_ELT(out, 3, Rf_ScalarString(STRING_ELT(g.lev, tj - 1)));
  static const char *const on[] = {"k", "seg", "i", "j"};
  clab::nlist(out, on);
  UNPROTECT(7);
  return klass(out, "crossk");
}

SEXP C_clab_r_qdrtst(SEXP x, SEXP y, SEXP win, SEXP nx, SEXP ny) {
  clab::who = "qdrtst";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  double wn[4];
  csp::win(win, wn);
  for (int a = 0; a < n; ++a) {
    if (!csp::inside(wn, p.x[a], p.y[a])) fail("every point must lie inside `win`.");
  }
  const int ax = clab::whole(nx, "nx", 1), ay = clab::whole(ny, "ny", 1);
  if (static_cast<double>(ax) * ay < 2 || static_cast<double>(ax) * ay > 1e7) {
    fail("`nx` times `ny` must be from 2 to 1e7.");
  }
  SEXP cn = PROTECT(Rf_allocMatrix(REALSXP, ay, ax));
  double st[3];
  kc(clab_qdrtst(p.x, p.y, n, wn, ax, ay, REAL(cn), st), "the kernel refused the input.");
  const double up = Rf_pchisq(st[0], st[1], 0, 0);
  static const char *const sn[] = {"x2", "df", "p", "ptwo", "lambda"};
  const double v[5] = {st[0], st[1], up, std::fmin(1.0, 2 * std::fmin(up, 1 - up)), st[2]};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(out, 0, cn);
  SET_VECTOR_ELT(out, 1, clab::nvec(5, sn, v));
  static const char *const on[] = {"counts", "stat"};
  clab::nlist(out, on);
  UNPROTECT(3);
  return klass(out, "qdrtst");
}

SEXP C_clab_r_morsta(SEXP counts, SEXP crit) {
  clab::who = "morsta";
  SEXP cv = mnum(counts, "counts", 0.0);
  const int n = clab::ilen(static_cast<double>(XLENGTH(cv)), "counts");
  if (n < 2) fail("at least 2 quadrats are needed.");
  const double a = clab::one(crit, "crit", 0.0, 1.0);
  double o[5];
  kc(clab_morsta(REAL(cv), n, Rf_qchisq(1 - a / 2, n - 1, 0, 0), Rf_qchisq(a / 2, n - 1, 0, 0), o),
     "the counts must add up to more than 1.");
  static const char *const nm[] = {"imor", "mclu", "muni", "imst", "chi", "p"};
  const double v[6] = {o[0], o[1], o[2], o[3], o[4], Rf_pchisq(o[4], n - 1, 0, 0)};
  UNPROTECT(1);
  return clab::nvec(6, nm, v);
}

SEXP C_clab_r_dissix(SEXP x, SEXP delta, SEXP w) {
  clab::who = "dissix";
  if (!Rf_isMatrix(x)) fail("`x` must be a matrix or data frame of counts, units by groups.");
  const int n = Rf_nrows(x), g = Rf_ncols(x);
  if (g < 2) fail("`x` must have at least 2 groups (columns).");
  SEXP xv = mnum(x, "x", 0.0);
  const double dl = clab::one(delta, "delta", 0.0, 1.0);
  SEXP gn = PROTECT(cnames(x, g, "g"));
  const bool hw = !Rf_isNull(w);
  wcsr c = hw ? wread(w, n) : wcsr{PROTECT(R_NilValue), nullptr, nullptr, nullptr};
  SEXP d = PROTECT(Rf_allocMatrix(REALSXP, g, g));
  SEXP dm = PROTECT(hw ? Rf_allocMatrix(REALSXP, g, g) : R_NilValue);
  SEXP gr = PROTECT(real(4 * static_cast<R_xlen_t>(g)));
  kc(clab_dissix(REAL(xv), n, g, dl, c.rp, c.ci, c.w, REAL(d), REAL(gr), hw ? REAL(dm) : nullptr),
     "every group needs a positive total, and `w` a positive weight.");
  SEXP dn = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dn, 0, gn);
  SET_VECTOR_ELT(dn, 1, gn);
  Rf_setAttrib(d, R_DimNamesSymbol, dn);
  if (hw) Rf_setAttrib(dm, R_DimNamesSymbol, dn);
  static const char *const cn[] = {"group", "is", "h", "atk", "gini"};
  SEXP cols = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(cols, 0, gn);
  for (int j = 0; j < 4; ++j) {
    SEXP v = real(g);
    SET_VECTOR_ELT(cols, j + 1, v);
    for (int k = 0; k < g; ++k) REAL(v)[k] = REAL(gr)[4 * static_cast<R_xlen_t>(k) + j];
  }
  SEXP gf = PROTECT(clab::frame(cols, cn));
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, d);
  SET_VECTOR_ELT(out, 1, gf);
  SET_VECTOR_ELT(out, 2, dm);
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(dl));
  static const char *const on[] = {"d", "grp", "dm", "delta"};
  clab::nlist(out, on);
  UNPROTECT(10);
  return klass(out, "dissix");
}

// ------------------------------------------------- geographic utilities

SEXP C_clab_r_pinpl(SEXP x, SEXP y, SEXP poly, SEXP pop) {
  clab::who = "pinpl";
  csp::pts p = csp::xy(x, y, R_NilValue, csp::planar);
  const int n = p.n;
  if (TYPEOF(poly) != VECSXP) fail("`poly` must be a data frame of x, y and id (and ring).");
  SEXP px = elt(poly, "x"), py = elt(poly, "y");
  if (Rf_isNull(px)) px = elt(poly, "lon");
  if (Rf_isNull(py)) py = elt(poly, "lat");
  SEXP pi = elt(poly, "id"), pr = elt(poly, "ring");
  if (Rf_isNull(px) || Rf_isNull(py) || Rf_isNull(pi)) {
    fail("`poly` must be a data frame of x, y and id (and ring).");
  }
  SEXP vx = PROTECT(num(px, "poly$x"));
  const int nv = clab::ilen(static_cast<double>(XLENGTH(vx)), "poly$x");
  if (nv < 3) fail("`poly` needs at least 3 vertices.");
  SEXP vy = PROTECT(num(py, "poly$y", nv));
  csp::lab id = csp::codes(pi, "poly$id", nv);
  SEXP rr;
  if (Rf_isNull(pr)) {
    rr = PROTECT(Rf_allocVector(INTSXP, nv));
    std::fill(INTEGER(rr), INTEGER(rr) + nv, 0);
    PROTECT(R_NilValue);
  } else {
    rr = csp::codes(pr, "poly$ring", nv).code;
  }
  // vertices grouped by polygon then ring, in their order within a ring
  std::vector<int> o(nv);
  for (int v = 0; v < nv; ++v) o[v] = v;
  const int *ic = INTEGER(id.code), *rc = INTEGER(rr);
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) {
    return ic[a] < ic[b] || (ic[a] == ic[b] && rc[a] < rc[b]);
  });
  std::vector<double> sx(nv), sy(nv);
  std::vector<int> sp(nv), sr(nv);
  for (int v = 0; v < nv; ++v) {
    sx[v] = REAL(vx)[o[v]];
    sy[v] = REAL(vy)[o[v]];
    sp[v] = ic[o[v]];
    sr[v] = rc[o[v]];
  }
  const int np = id.g;
  SEXP pv = PROTECT(Rf_isNull(pop) ? R_NilValue : num(pop, "pop", np, 0.0));
  SEXP reg = PROTECT(Rf_allocVector(INTSXP, n));
  SEXP nreg = PROTECT(Rf_allocVector(INTSXP, n));
  kc(clab_pinpl(p.x, p.y, n, sx.data(), sy.data(), sp.data(), sr.data(), nv, np, INTEGER(reg),
                INTEGER(nreg)),
     "the kernel refused the input.");
  SEXP lab = PROTECT(Rf_allocVector(STRSXP, n));
  std::vector<double> units(np, 0.0);
  for (int i = 0; i < n; ++i) {
    const int g = INTEGER(reg)[i];
    SET_STRING_ELT(lab, i, g >= 0 ? STRING_ELT(id.lev, g) : NA_STRING);
    if (g >= 0) units[g] += 1;
  }
  // the regions, and the population of those holding a unit
  double tot = 0, cov = 0, with = 0, nu = 0;
  const double na = std::numeric_limits<double>::quiet_NaN();
  for (int q = 0; q < np; ++q) {
    const double pp = Rf_isNull(pv) ? na : REAL(pv)[q];
    tot += pp;
    if (units[q] > 0) {
      cov += pp;
      with += 1;
    }
    nu += units[q];
  }
  static const char *const rn[] = {"region", "units", "pop", "hasunit", "popshr"};
  SEXP cols = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(cols, 0, id.lev);
  SEXP cu = real(np);
  SET_VECTOR_ELT(cols, 1, cu);
  SEXP cpop = real(np);
  SET_VECTOR_ELT(cols, 2, cpop);
  SEXP ch = Rf_allocVector(LGLSXP, np);
  SET_VECTOR_ELT(cols, 3, ch);
  SEXP cs = real(np);
  SET_VECTOR_ELT(cols, 4, cs);
  for (int q = 0; q < np; ++q) {
    REAL(cu)[q] = units[q];
    REAL(cpop)[q] = Rf_isNull(pv) ? na : REAL(pv)[q];
    LOGICAL(ch)[q] = units[q] > 0;
    REAL(cs)[q] = tot > 0 ? REAL(cpop)[q] / tot : na;
  }
  SEXP rf = PROTECT(clab::frame(cols, rn));
  static const char *const cn[] = {"regions", "withunit", "units", "totpop",
                                   "covpop",  "uncpop",   "covshr"};
  const double cv[7] = {static_cast<double>(np), with, nu, tot, cov, tot - cov,
                        tot > 0 ? cov / tot : na};
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 4));
  SET_VECTOR_ELT(out, 0, lab);
  SET_VECTOR_ELT(out, 1, nreg);
  SET_VECTOR_ELT(out, 2, rf);
  SET_VECTOR_ELT(out, 3, clab::nvec(7, cn, cv));
  static const char *const on[] = {"reg", "nreg", "regs", "cov"};
  clab::nlist(out, on);
  UNPROTECT(14);
  return klass(out, "pinpl");
}

SEXP C_clab_r_arwt(SEXP x, SEXP id, SEXP from, SEXP to, SEXP area, SEXP top) {
  clab::who = "arwt";
  SEXP xv = PROTECT(num(x, "x"));
  const int ns = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  csp::lab src = csp::codes(id, "id", ns);
  if (src.g != ns) fail("`id` must name each source unit once.");
  const int np = clab::ilen(static_cast<double>(Rf_xlength(from)), "from");
  if (np < 1) fail("`from` must hold at least one pair.");
  csp::lab tg = csp::codes(to, "to", np);
  csp::lab fr = csp::codes(from, "from", np);
  // the pairs' sources as positions in id
  std::vector<int> fi(np);
  for (int q = 0; q < np; ++q) {
    SEXP lb = STRING_ELT(fr.lev, INTEGER(fr.code)[q]);
    int k = -1;
    for (int s = 0; s < src.g && k < 0; ++s) {
      if (clab::seql(STRING_ELT(src.lev, s), lb)) k = s;
    }
    if (k < 0) fail("`from` names a unit not in `id`: %s.", Rf_translateCharUTF8(lb));
    fi[q] = k;
  }
  // the values in the order of the sorted labels
  std::vector<double> val(ns);
  for (int s = 0; s < ns; ++s) val[INTEGER(src.code)[s]] = REAL(xv)[s];
  SEXP av = PROTECT(Rf_isNull(area) ? R_NilValue : num(area, "area", np, 0.0));
  std::vector<double> ar(np, 1.0);
  if (!Rf_isNull(av)) std::copy(REAL(av), REAL(av) + np, ar.begin());
  SEXP tv = PROTECT(num(top, "top", -1, 1.0));
  const int nt = tg.g;
  SEXP ov = PROTECT(real(nt));
  std::vector<double> sh(np);
  std::vector<int> mn(ns);
  double lost = 0;
  kc(clab_arwt(val.data(), ns, fi.data(), INTEGER(tg.code), ar.data(), np, nt, REAL(ov), sh.data(),
               mn.data(), &lost),
     "the kernel refused the input.");
  SEXP oc = PROTECT(Rf_allocVector(VECSXP, 2));
  SET_VECTOR_ELT(oc, 0, tg.lev);
  SET_VECTOR_ELT(oc, 1, ov);
  static const char *const onm[] = {"to", "value"};
  SEXP of = PROTECT(clab::frame(oc, onm));
  // each source's main target and its share
  SEXP mc = PROTECT(Rf_allocVector(VECSXP, 3));
  SET_VECTOR_ELT(mc, 0, src.lev);
  SEXP mt = Rf_allocVector(STRSXP, ns);
  SET_VECTOR_ELT(mc, 1, mt);
  SEXP ms = real(ns);
  SET_VECTOR_ELT(mc, 2, ms);
  std::vector<double> best(ns, -1.0);
  for (int q = 0; q < np; ++q) best[fi[q]] = std::max(best[fi[q]], sh[q]);
  bool exact = true;
  std::vector<int> cnt(ns, 0);
  for (int q = 0; q < np; ++q) ++cnt[fi[q]];
  for (int s = 0; s < ns; ++s) {
    SET_STRING_ELT(mt, s, mn[s] >= 0 ? STRING_ELT(tg.lev, mn[s]) : NA_STRING);
    REAL(ms)[s] = mn[s] >= 0 ? best[s] : std::numeric_limits<double>::quiet_NaN();
    exact = exact && cnt[s] <= 1;
  }
  static const char *const mnm[] = {"from", "to", "share"};
  SEXP mf = PROTECT(clab::frame(mc, mnm));
  // concentration of the targets' values: Gini and the shares of the top
  // units
  std::vector<double> v(REAL(ov), REAL(ov) + nt);
  std::sort(v.begin(), v.end());
  double sum = 0, cum = 0, scum = 0;
  for (int q = 0; q < nt; ++q) {
    sum += v[q];
  }
  for (int q = 0; q < nt; ++q) {
    cum += v[q];
    scum += cum;
  }
  const R_xlen_t ntop = XLENGTH(tv);
  std::vector<double> cc(1 + ntop);
  cc[0] = sum > 0 ? (nt + 1 - 2 * scum / sum) / nt : 0.0;
  for (R_xlen_t k = 0; k < ntop; ++k) {
    const int m = static_cast<int>(std::min<double>(nt, std::floor(REAL(tv)[k])));
    double s = 0;
    for (int q = nt - m; q < nt; ++q) s += v[q];
    cc[1 + k] = sum > 0 ? s / sum : std::numeric_limits<double>::quiet_NaN();
  }
  SEXP conc = PROTECT(real(1 + ntop));
  SEXP cnn = PROTECT(Rf_allocVector(STRSXP, 1 + ntop));
  SET_STRING_ELT(cnn, 0, Rf_mkChar("gini"));
  char buf[48];
  for (R_xlen_t k = 0; k <= ntop; ++k) {
    REAL(conc)[k] = cc[k];
    if (k > 0) {
      std::snprintf(buf, sizeof buf, "top%g", std::floor(REAL(tv)[k - 1]));
      SET_STRING_ELT(cnn, k, Rf_mkChar(buf));
    }
  }
  Rf_setAttrib(conc, R_NamesSymbol, cnn);
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, of);
  SET_VECTOR_ELT(out, 1, mf);
  SET_VECTOR_ELT(out, 2, Rf_ScalarLogical(exact));
  SET_VECTOR_ELT(out, 3, Rf_ScalarReal(lost));
  SET_VECTOR_ELT(out, 4, conc);
  static const char *const on[] = {"out", "map", "exact", "lost", "conc"};
  clab::nlist(out, on);
  UNPROTECT(17);
  return klass(out, "arwt");
}

}  // extern "C"
