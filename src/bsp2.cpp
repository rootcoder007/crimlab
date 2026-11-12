// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The R entry points of the spatial weights and the statistics built on
// them (Geary's C, the global and local G of Getis and Ord, local Moran)
// and of the directional mean. Weights come from spwt() or as an spdep
// "listw" object; both are turned into compressed sparse rows here.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <vector>

#include "crimlab.h"
#include "rchk.h"
#include "spin.h"

using clab::fail;
using clab::num;

namespace {

SEXP real(R_xlen_t n) { return Rf_allocVector(REALSXP, n); }
SEXP ints(R_xlen_t n) { return Rf_allocVector(INTSXP, n); }

void kc(int rc, const char *msg) {
  if (rc != CLAB_OK) fail("%s", msg);
}

// an element of list x by name, or R_NilValue
SEXP elt(SEXP x, const char *nm) {
  SEXP n = Rf_getAttrib(x, R_NamesSymbol);
  for (R_xlen_t j = 0; j < Rf_xlength(n); ++j) {
    if (std::strcmp(CHAR(STRING_ELT(n, j)), nm) == 0) return VECTOR_ELT(x, j);
  }
  return R_NilValue;
}

// Weights as compressed sparse rows: rp (n + 1), ci and w (nnz), rows in
// increasing order of neighbour. v holds the three vectors, protected by
// the caller.
struct csr {
  SEXP v;
  int n;
  const int *rp, *ci;
  const double *w;
};

// edges (i, j, w), 0-based, into rows sorted by j; duplicates and self
// edges refused
csr build(int n, const std::vector<int> &ei, const std::vector<int> &ej,
          const std::vector<double> &ew) {
  const std::size_t m = ei.size();
  std::vector<int> o(m);
  for (std::size_t q = 0; q < m; ++q) o[q] = static_cast<int>(q);
  std::sort(o.begin(), o.end(), [&](int a, int b) {
    return ei[a] < ei[b] || (ei[a] == ei[b] && ej[a] < ej[b]);
  });
  csr c{PROTECT(Rf_allocVector(VECSXP, 3)), n, nullptr, nullptr, nullptr};
  SET_VECTOR_ELT(c.v, 0, ints(n + 1));
  SET_VECTOR_ELT(c.v, 1, ints(static_cast<R_xlen_t>(m)));
  SET_VECTOR_ELT(c.v, 2, real(static_cast<R_xlen_t>(m)));
  int *rp = INTEGER(VECTOR_ELT(c.v, 0)), *ci = INTEGER(VECTOR_ELT(c.v, 1));
  double *w = REAL(VECTOR_ELT(c.v, 2));
  std::fill(rp, rp + n + 1, 0);
  for (std::size_t q = 0; q < m; ++q) {
    const int a = o[q];
    if (ei[a] == ej[a]) fail("a unit cannot be its own neighbour in `w`.");
    if (q > 0 && ei[a] == ei[o[q - 1]] && ej[a] == ej[o[q - 1]]) {
      fail("`w` lists a pair of units twice.");
    }
    ci[q] = ej[a];
    w[q] = ew[a];
    ++rp[ei[a] + 1];
  }
  for (int i = 0; i < n; ++i) rp[i + 1] += rp[i];
  c.rp = rp;
  c.ci = ci;
  c.w = w;
  return c;
}

// the weights argument: spwt()'s list or an spdep listw, for n units
csr wread(SEXP w, int n) {
  std::vector<int> ei, ej;
  std::vector<double> ew;
  if (Rf_inherits(w, "listw")) {
    SEXP nb = elt(w, "neighbours"), wt = elt(w, "weights");
    if (TYPEOF(nb) != VECSXP || TYPEOF(wt) != VECSXP || XLENGTH(nb) != n ||
        XLENGTH(wt) != n) {
      fail("`w` must have one row of neighbours and weights per value of `x`.");
    }
    for (int i = 0; i < n; ++i) {
      SEXP a = VECTOR_ELT(nb, i), b = VECTOR_ELT(wt, i);
      if (TYPEOF(a) != INTSXP) fail("`w` must have integer neighbour lists.");
      const R_xlen_t k = XLENGTH(a);
      if (k == 1 && INTEGER_RO(a)[0] == 0) continue;  // no neighbours
      if (TYPEOF(b) != REALSXP || XLENGTH(b) != k) {
        fail("`w` must have one weight per neighbour.");
      }
      for (R_xlen_t q = 0; q < k; ++q) {
        const int j = INTEGER_RO(a)[q];
        if (j < 1 || j > n) fail("`w` names a neighbour outside 1 to %d.", n);
        ei.push_back(i);
        ej.push_back(j - 1);
        ew.push_back(REAL_RO(b)[q]);
      }
    }
  } else {
    SEXP nn = TYPEOF(w) == VECSXP ? elt(w, "n") : R_NilValue;
    if (Rf_isNull(nn)) fail("`w` must be the list spwt() returns or an spdep listw.");
    if (clab::val1(nn) != n) fail("`w` is for %g units, not %d.", clab::val1(nn), n);
    SEXP i = PROTECT(num(elt(w, "i"), "w$i", -1, 1.0));
    const R_xlen_t m = XLENGTH(i);
    SEXP j = PROTECT(num(elt(w, "j"), "w$j", m, 1.0));
    SEXP v = PROTECT(num(elt(w, "w"), "w$w", m, 0.0));
    for (R_xlen_t q = 0; q < m; ++q) {
      const double a = REAL_RO(i)[q], b = REAL_RO(j)[q];
      if (a > n || b > n || a != std::floor(a) || b != std::floor(b)) {
        fail("`w$i` and `w$j` must be unit numbers from 1 to %d.", n);
      }
      ei.push_back(static_cast<int>(a) - 1);
      ej.push_back(static_cast<int>(b) - 1);
      ew.push_back(REAL_RO(v)[q]);
    }
    UNPROTECT(3);
  }
  return build(n, ei, ej, ew);
}

}  // namespace

extern "C" {

// --------------------------------------------------------------- weights

SEXP C_clab_r_spwt(SEXP x, SEXP y, SEXP type, SEXP k, SEXP d, SEXP nb,
                   SEXP idw, SEXP style, SEXP metric, SEXP lonlat) {
  clab::who = "spwt";
  static const char *const to[] = {"knn", "band", "contig"};
  const int tp = csp::pick(type, "type", to, 3);
  static const char *const so[] = {"row", "raw"};
  const int rs = csp::pick(style, "style", so, 2) == 0;
  static const char *const mo[] = {"euclid", "manhattan"};
  const int me = csp::pick(metric, "metric", mo, 2);
  const double pw = clab::val1(idw);
  if (!(pw >= 0) || !std::isfinite(pw)) fail("`idw` must be one number from 0.");
  const bool pts = !Rf_isNull(x);
  csp::pts p{};
  int np = 0;
  if (pts) {
    p = csp::xy(x, y, lonlat);
    np = 1;
  }
  int n;
  SEXP rpv, civ, dv;
  if (tp == 2) {
    if (TYPEOF(nb) != VECSXP) fail("`nb` must be a list of neighbour numbers, one element per unit.");
    n = clab::ilen(static_cast<double>(XLENGTH(nb)), "nb");
    if (pts && p.n != n) fail("`nb` must have one element per point.");
    if (!pts && pw > 0) fail("inverse-distance weights need the points `x`.");
    std::vector<int> ei, ej;
    for (int i = 0; i < n; ++i) {
      SEXP a = VECTOR_ELT(nb, i);
      if (!clab::isnum(a)) fail("`nb` must be a list of neighbour numbers, one element per unit.");
      for (R_xlen_t q = 0; q < XLENGTH(a); ++q) {
        const double j = TYPEOF(a) == INTSXP ? INTEGER_RO(a)[q] : REAL_RO(a)[q];
        if (j == 0 && XLENGTH(a) == 1) break;  // spdep's "no neighbours"
        if (!(j >= 1 && j <= n && j == std::floor(j))) {
          fail("`nb` must hold unit numbers from 1 to %d.", n);
        }
        ei.push_back(i);
        ej.push_back(static_cast<int>(j) - 1);
      }
    }
    std::vector<double> ed(ei.size());
    for (std::size_t q = 0; q < ei.size(); ++q) {
      ed[q] = pts ? (me ? std::fabs(p.x[ei[q]] - p.x[ej[q]]) + std::fabs(p.y[ei[q]] - p.y[ej[q]])
                        : std::hypot(p.x[ei[q]] - p.x[ej[q]], p.y[ei[q]] - p.y[ej[q]]))
                  : clab::ninf;
    }
    csr c = build(n, ei, ej, ed);
    rpv = VECTOR_ELT(c.v, 0);
    civ = VECTOR_ELT(c.v, 1);
    dv = VECTOR_ELT(c.v, 2);
    if (!pts) {
      for (R_xlen_t q = 0; q < XLENGTH(dv); ++q) REAL(dv)[q] = NA_REAL;
    }
    np += 1;
  } else {
    if (!pts) fail("`x` must give the points for `type = \"%s\"`.", to[tp]);
    n = p.n;
    if (n < 2) fail("at least 2 points are needed.");
    SEXP iw = PROTECT(ints(2 * static_cast<R_xlen_t>(n) + 2));
    if (tp == 0) {
      const int kk = clab::whole(k, "k", 1);
      if (kk > n - 1) fail("`k` must be below the number of points.");
      const R_xlen_t m = static_cast<R_xlen_t>(n) * kk;
      rpv = PROTECT(ints(n + 1));
      civ = PROTECT(ints(m));
      dv = PROTECT(real(m));
      kc(clab_swknb(p.x, p.y, n, kk, me, INTEGER(civ), REAL(dv), INTEGER(iw)),
         "the kernel refused the input.");
      for (int i = 0; i <= n; ++i) INTEGER(rpv)[i] = i * kk;
      // rows in increasing order of neighbour, as the statistics need
      for (int i = 0; i < n; ++i) {
        int *c = INTEGER(civ) + static_cast<R_xlen_t>(i) * kk;
        double *e = REAL(dv) + static_cast<R_xlen_t>(i) * kk;
        for (int u = 1; u < kk; ++u) {
          const int j = c[u];
          const double f = e[u];
          int v = u - 1;
          while (v >= 0 && c[v] > j) {
            c[v + 1] = c[v];
            e[v + 1] = e[v];
            --v;
          }
          c[v + 1] = j;
          e[v + 1] = f;
        }
      }
    } else {
      const double dm = clab::one(d, "d", 0.0, std::numeric_limits<double>::infinity());
      int nnz = 0;
      rpv = PROTECT(ints(n + 1));
      kc(clab_swband(p.x, p.y, n, dm, me, 0, INTEGER(rpv), nullptr, nullptr, &nnz,
                     INTEGER(iw)),
         "the kernel refused the input.");
      civ = PROTECT(ints(nnz));
      dv = PROTECT(real(nnz));
      kc(clab_swband(p.x, p.y, n, dm, me, nnz, INTEGER(rpv), INTEGER(civ), REAL(dv),
                     &nnz, INTEGER(iw)),
         "the kernel refused the input.");
    }
    np += 4;
  }
  const R_xlen_t m = XLENGTH(civ);
  SEXP wv = PROTECT(real(m));
  kc(clab_swwt(INTEGER(rpv), REAL(dv), n, pw, rs, REAL(wv)),
     "inverse-distance weights need distinct points: two neighbours coincide.");
  SEXP iv = PROTECT(ints(m));
  SEXP jv = PROTECT(ints(m));
  for (int i = 0; i < n; ++i) {
    for (int q = INTEGER(rpv)[i]; q < INTEGER(rpv)[i + 1]; ++q) {
      INTEGER(iv)[q] = i + 1;
      INTEGER(jv)[q] = INTEGER(civ)[q] + 1;
    }
  }
  SEXP out = PROTECT(Rf_allocVector(VECSXP, 5));
  SET_VECTOR_ELT(out, 0, Rf_ScalarInteger(n));
  SET_VECTOR_ELT(out, 1, iv);
  SET_VECTOR_ELT(out, 2, jv);
  SET_VECTOR_ELT(out, 3, dv);
  SET_VECTOR_ELT(out, 4, wv);
  static const char *const nm[] = {"n", "i", "j", "d", "w"};
  UNPROTECT(np + 4);
  return clab::nlist(out, nm);
}

// ----------------------------------------------- spatial autocorrelation

SEXP C_clab_r_gearc(SEXP x, SEXP w, SEXP nsim, SEXP seed) {
  clab::who = "gearc";
  SEXP xv = PROTECT(num(x, "x"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  if (n < 4) fail("at least 4 values are needed.");
  csr c = wread(w, n);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int s = ns > 0 ? csp::seed(seed) : 0u;
  SEXP wk = PROTECT(real(n));
  double o[9];
  kc(clab_gearc(REAL_RO(xv), n, c.rp, c.ci, c.w, ns, s, REAL(wk), o),
     "`x` is constant, or `w` has no positive weight: C is undefined.");
  static const char *const nm[] = {"c", "vn", "vr", "zn", "zr", "pn", "pr", "pperm", "s0"};
  UNPROTECT(3);
  return clab::nvec(9, nm, o);
}

SEXP C_clab_r_gorg(SEXP x, SEXP w, SEXP nsim, SEXP seed) {
  clab::who = "gorg";
  SEXP xv = PROTECT(num(x, "x", -1, 0.0));
  const int n = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  if (n < 4) fail("at least 4 values are needed.");
  csr c = wread(w, n);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int s = ns > 0 ? csp::seed(seed) : 0u;
  SEXP wk = PROTECT(real(n));
  double o[6];
  kc(clab_gorg(REAL_RO(xv), n, c.rp, c.ci, c.w, ns, s, REAL(wk), o),
     "G is undefined: fewer than 2 positive values of `x`, or no positive weight.");
  static const char *const nm[] = {"g", "e", "v", "z", "p", "pperm"};
  UNPROTECT(3);
  return clab::nvec(6, nm, o);
}

SEXP C_clab_r_lmorn(SEXP x, SEXP w, SEXP nsim, SEXP seed) {
  clab::who = "lmorn";
  SEXP xv = PROTECT(num(x, "x"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  if (n < 3) fail("at least 3 values are needed.");
  csr c = wread(w, n);
  const int ns = clab::whole(nsim, "nsim", 0);
  const unsigned int s = ns > 0 ? csp::seed(seed) : 0u;
  SEXP iw = PROTECT(ints(n));
  SEXP wk = PROTECT(real(n));
  static const char *const nm[] = {"ii", "e", "v", "z", "p", "pperm", "quad"};
  double *col[7];
  SEXP v = PROTECT(Rf_allocVector(VECSXP, 7));
  for (int j = 0; j < 7; ++j) {
    SET_VECTOR_ELT(v, j, real(n));
    col[j] = REAL(VECTOR_ELT(v, j));
  }
  SEXP o = PROTECT(real(7 * static_cast<R_xlen_t>(n)));
  kc(clab_lmorn(REAL_RO(xv), n, c.rp, c.ci, c.w, ns, s, INTEGER(iw), REAL(wk), REAL(o)),
     "`x` is constant: local Moran is undefined.");
  for (int j = 0; j < 7; ++j) {
    std::memcpy(col[j], REAL(o) + static_cast<R_xlen_t>(j) * n, n * sizeof(double));
  }
  SEXP f = clab::frame(v, nm);
  UNPROTECT(6);
  return f;
}

SEXP C_clab_r_locgi(SEXP x, SEXP w, SEXP star) {
  clab::who = "locgi";
  SEXP xv = PROTECT(num(x, "x"));
  const int n = clab::ilen(static_cast<double>(XLENGTH(xv)), "x");
  if (n < 3) fail("at least 3 values are needed.");
  csr c = wread(w, n);
  const bool st = csp::flag(star, "star");
  SEXP out = PROTECT(real(n));
  kc(clab_locgi(REAL_RO(xv), n, c.rp, c.ci, c.w, st ? 1 : 0, REAL(out)),
     "the kernel refused the input.");
  UNPROTECT(3);
  return out;
}

// --------------------------------------------------------- directions

SEXP C_clab_r_dirmn(SEXP x, SEXP y, SEXP origin, SEXP w, SEXP angle,
                    SEXP lonlat) {
  clab::who = "dirmn";
  if (Rf_isNull(angle) == Rf_isNull(x)) {
    fail("give the points `x` (with `origin`) or the angles `angle`, not both.");
  }
  SEXP a;
  int n;
  int np = 0;
  if (!Rf_isNull(angle)) {
    a = PROTECT(clab::dcopy(PROTECT(num(angle, "angle"))));
    n = clab::ilen(static_cast<double>(XLENGTH(a)), "angle");
    np = 2;
  } else {
    csp::pts p = csp::xy(x, y, lonlat);
    SEXP o = PROTECT(num(origin, "origin", 2));
    double ox = REAL_RO(o)[0], oy = REAL_RO(o)[1];
    csp::fwd(p, &ox, &oy, 1);
    n = p.n;
    a = PROTECT(real(n));
    for (int i = 0; i < n; ++i) {
      const double dx = p.x[i] - ox, dy = p.y[i] - oy;
      if (dx == 0 && dy == 0) fail("a point lies on `origin`: it has no direction.");
      REAL(a)[i] = std::atan2(dy, dx);
    }
    np = 3;
  }
  const double *wp = nullptr;
  if (!Rf_isNull(w)) {
    SEXP wv = PROTECT(num(w, "w", n, 0.0));
    wp = REAL_RO(wv);
    ++np;
  }
  double o[5];
  kc(clab_dirmn(REAL(a), wp, n, o), "`w` must have a positive sum.");
  // the compass bearing (degrees clockwise from north) of the mean
  const double pi = 3.14159265358979323846;
  double b = 90.0 - o[0] * 180.0 / pi;
  if (b < 0) b += 360.0;
  const double r[7] = {o[0], b, o[1], o[2], o[3], o[4], static_cast<double>(n)};
  static const char *const nm[] = {"theta", "bearing", "r", "var", "sd", "w", "n"};
  UNPROTECT(np);
  return clab::nvec(7, nm, r);
}

}  // extern "C"
