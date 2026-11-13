// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Helpers of the spatial R entry points; see spin.h.
#include "spin.h"

#include <R_ext/Random.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "rchk.h"

namespace csp {

namespace {

constexpr double earth = 6371008.8;  // mean Earth radius, metres
constexpr double deg = 3.14159265358979323846 / 180.0;

// the position of the column named a in names nm, or -1
int find(SEXP nm, const char *a) {
  const R_xlen_t k = Rf_xlength(nm);
  for (R_xlen_t j = 0; j < k; ++j) {
    SEXP s = STRING_ELT(nm, j);
    if (s != NA_STRING && std::strcmp(CHAR(s), a) == 0) return static_cast<int>(j);
  }
  return -1;
}

// column j of the numeric matrix m with nr rows, as a fresh vector
SEXP mcol(SEXP m, R_xlen_t nr, int j) {
  SEXP v = PROTECT(Rf_allocVector(REALSXP, nr));
  double *o = REAL(v);
  const R_xlen_t at = nr * j;
  for (R_xlen_t i = 0; i < nr; ++i) {
    if (TYPEOF(m) == REALSXP) {
      o[i] = REAL_RO(m)[at + i];
    } else {
      const int k = INTEGER_RO(m)[at + i];
      o[i] = k == NA_INTEGER ? NA_REAL : k;
    }
  }
  UNPROTECT(1);
  return v;
}

// The kind of an sf coordinate reference system: 1 geographic, 0
// projected, -1 unknown (none, or no WKT). WKT2 names a geographic CRS
// GEOGCRS (or GEODCRS), WKT1 GEOGCS.
int crskind(SEXP crs) {
  if (TYPEOF(crs) != VECSXP) return -1;
  const int j = find(Rf_getAttrib(crs, R_NamesSymbol), "wkt");
  SEXP w = j < 0 ? R_NilValue : VECTOR_ELT(crs, j);
  if (TYPEOF(w) != STRSXP || XLENGTH(w) != 1 || STRING_ELT(w, 0) == NA_STRING) {
    return -1;
  }
  const char *s = CHAR(STRING_ELT(w, 0));
  while (*s == ' ' || *s == '\n') ++s;
  return std::strncmp(s, "GEOGCRS", 7) == 0 || std::strncmp(s, "GEODCRS", 7) == 0 ||
         std::strncmp(s, "GEOGCS", 6) == 0;
}

}  // namespace

//' @srrstats {G2.7} xy() takes coordinates as two vectors, a matrix or a
//'   data frame (base or a subclass such as a tibble, read as its list of
//'   columns).
//' @srrstats {G2.8} xy() is the one pre-processing routine of every
//'   spatial function: it turns each accepted form into two double vectors
//'   in planar units before any kernel sees them.
//' @srrstats {G2.10} A coordinate column is found by its name, the same
//'   way for a matrix, a base data frame and a tibble; no column is taken
//'   by position when names exist.
//' @srrstats {G2.11} A column with a class of its own (units, say) is read
//'   as its numbers, as other classed numbers are (G2.6); test-sp.R.
//' @srrstats {G2.12} A list column is refused with an error naming it
//'   (test-sp.R).
//' @srrstats {SP2.7} xy() validates the coordinates' form, names, values
//'   and, for longitude and latitude, range.
//' @srrstats {SP2.8} xy() is the single pre-processing routine; every
//'   spatial entry point calls it first.
//' @srrstats {SP2.0b} Longitude and latitude are never guessed: columns
//'   named lon and lat need `lonlat = TRUE`, longitudes outside -180 to 180
//'   and latitudes outside -90 to 90 stop with an error, functions that
//'   take a rectangle in planar units refuse lon/lat columns, and an sf
//'   object's coordinate reference system must agree with `lonlat`.
//' @srrstats {SP2.0} Points come as an sf or sfc object of POINT
//'   geometries (the class of the sf package) or as coordinate vectors,
//'   matrices and data frames (the form sf::st_coordinates() returns).
//' @srrstats {SP2.1} The sp package is not used; sf objects are read
//'   directly.
//' @srrstats {SP2.4} An sf object's coordinate reference system is read
//'   from its WKT2 representation (the `wkt` field), never from a PROJ
//'   string.
//' @srrstats {SP2.4a} No function takes a PROJ string; the CRS of sf input
//'   is read as WKT.
//' @noRd
pts xy(SEXP x, SEXP y, SEXP lonlat, mode m) {
  pts p{};
  p.ll = m == raw || (m == proj && flag(lonlat, "lonlat"));
  SEXP cx = x, cy = y;
  const char *nx = p.ll ? "lon" : "x", *ny = p.ll ? "lat" : "y";
  int np = 0;
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (Rf_isNull(y) && m == raw && clab::isnum(x) && Rf_isNull(d) &&
      XLENGTH(x) == 2) {
    // one point, c(lon, lat)
    cx = PROTECT(Rf_ScalarReal(Rf_asReal(x)));
    cy = PROTECT(Rf_ScalarReal(REAL_RO(PROTECT(clab::dbl(x)))[1]));
    np = 3;
  } else if (Rf_isNull(y) && (Rf_inherits(x, "sf") || Rf_inherits(x, "sfc"))) {
    SEXP g = x;
    if (Rf_inherits(x, "sf")) {
      SEXP sc = Rf_getAttrib(x, Rf_install("sf_column"));
      const int j = TYPEOF(sc) == STRSXP && XLENGTH(sc) == 1
                        ? find(Rf_getAttrib(x, R_NamesSymbol), CHAR(STRING_ELT(sc, 0)))
                        : -1;
      g = j < 0 ? R_NilValue : VECTOR_ELT(x, j);
    }
    if (!Rf_inherits(g, "sfc_POINT")) {
      clab::fail("`x` must hold POINT geometries.");
    }
    const int geo = crskind(Rf_getAttrib(g, Rf_install("crs")));
    if (geo == 1 && m == planar) {
      clab::fail("`x` has a geographic CRS: this function takes planar "
                 "coordinates; project it first (eqproj() or "
                 "sf::st_transform()).");
    }
    if (geo == 1 && !p.ll) clab::fail("`x` has a geographic CRS: set `lonlat = TRUE`.");
    if (geo == 0 && p.ll) {
      clab::fail("`x` has a projected CRS: its coordinates are not longitude "
                 "and latitude.");
    }
    const R_xlen_t n = XLENGTH(g);
    cx = PROTECT(Rf_allocVector(REALSXP, n));
    cy = PROTECT(Rf_allocVector(REALSXP, n));
    np = 2;
    for (R_xlen_t i = 0; i < n; ++i) {
      SEXP e = VECTOR_ELT(g, i);
      if (TYPEOF(e) != REALSXP || XLENGTH(e) != 2) {
        clab::fail("`x` must hold two-dimensional points.");
      }
      REAL(cx)[i] = REAL_RO(e)[0];
      REAL(cy)[i] = REAL_RO(e)[1];
    }
  } else if (Rf_isNull(y)) {
    SEXP nm = R_NilValue;
    R_xlen_t nc = 0, nr = 0;
    const bool lst = TYPEOF(x) == VECSXP;
    if (lst) {
      nm = Rf_getAttrib(x, R_NamesSymbol);
      nc = XLENGTH(x);
    } else if (clab::isnum(x) && TYPEOF(d) == INTSXP && XLENGTH(d) == 2) {
      nr = INTEGER_RO(d)[0];
      nc = INTEGER_RO(d)[1];
      SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
      if (!Rf_isNull(dn)) nm = VECTOR_ELT(dn, 1);
    } else {
      clab::fail("`y` is missing: give `x` and `y`, or a matrix or data "
                 "frame of coordinates as `x`.");
    }
    int jx = -1, jy = -1;
    bool isll = false;
    if (!Rf_isNull(nm)) {
      if (m != raw) {
        jx = find(nm, "x");
        jy = find(nm, "y");
        if (jx < 0 || jy < 0) {
          jx = find(nm, "X");
          jy = find(nm, "Y");
        }
      }
      if (jx < 0 || jy < 0) {
        jx = find(nm, "lon");
        jy = find(nm, "lat");
        isll = true;
      }
    } else if (nc == 2) {
      jx = 0;
      jy = 1;
    }
    if (jx < 0 || jy < 0) {
      clab::fail(m == raw ? "`x` must have columns named lon and lat, or be "
                            "a matrix of exactly two unnamed columns."
                          : "`x` must have columns named x and y, or lon and "
                            "lat, or be a matrix of exactly two unnamed "
                            "columns.");
    }
    if (isll && m == planar) {
      clab::fail("columns lon and lat hold degrees: this function takes "
                 "planar coordinates; project them first with eqproj().");
    }
    if (isll && !p.ll) {
      clab::fail("columns lon and lat hold degrees: set `lonlat = TRUE`.");
    }
    if (lst) {
      cx = VECTOR_ELT(x, jx);
      cy = VECTOR_ELT(x, jy);
    } else {
      cx = PROTECT(mcol(x, nr, jx));
      cy = PROTECT(mcol(x, nr, jy));
      np = 2;
    }
  }
  SEXP xyv[2] = {cx, cy};
  static const char *const xyn[] = {"x", "y"};
  const int nb = byname(xyv, xyn, 2);
  np += nb;
  cy = xyv[1];
  SEXP vx = PROTECT(clab::num(cx, nx));
  const int n = clab::ilen(static_cast<double>(XLENGTH(vx)), nx);
  SEXP vy = PROTECT(clab::num(cy, ny, n));
  SEXP v = PROTECT(Rf_allocVector(VECSXP, 2));
  if (!p.ll) {
    SET_VECTOR_ELT(v, 0, clab::dcopy(vx));
    SET_VECTOR_ELT(v, 1, clab::dcopy(vy));
  } else {
    const double *lo = REAL_RO(vx), *la = REAL_RO(vy);
    double ss = 0.0, sc = 0.0, sl = 0.0;
    for (int i = 0; i < n; ++i) {
      if (lo[i] < -180.0 || lo[i] > 180.0) {
        clab::fail("`%s` must be longitudes from -180 to 180.", nx);
      }
      if (la[i] < -90.0 || la[i] > 90.0) {
        clab::fail("`%s` must be latitudes from -90 to 90.", ny);
      }
      ss += std::sin(lo[i] * deg);
      sc += std::cos(lo[i] * deg);
      sl += la[i];
    }
    p.lon0 = std::atan2(ss, sc) / deg;
    p.lat0 = sl / n;
    if (m == raw) {
      SET_VECTOR_ELT(v, 0, clab::dcopy(vx));
      SET_VECTOR_ELT(v, 1, clab::dcopy(vy));
    } else {
      SET_VECTOR_ELT(v, 0, Rf_allocVector(REALSXP, n));
      SET_VECTOR_ELT(v, 1, Rf_allocVector(REALSXP, n));
      if (clab_eqproj(lo, la, n, p.lon0, p.lat0, earth, REAL(VECTOR_ELT(v, 0)),
                      REAL(VECTOR_ELT(v, 1))) != CLAB_OK) {
        clab::fail("the points reach the antipode of their centre: no "
                   "equal-area projection about it holds them all.");
      }
    }
  }
  UNPROTECT(3 + np);
  PROTECT(v);
  p.v = v;
  p.x = REAL_RO(VECTOR_ELT(v, 0));
  p.y = REAL_RO(VECTOR_ELT(v, 1));
  p.n = n;
  return p;
}

void fwd(const pts &p, double *px, double *py, int n) {
  if (!p.ll) return;
  for (int i = 0; i < n; ++i) {
    if (px[i] < -180.0 || px[i] > 180.0 || py[i] < -90.0 || py[i] > 90.0) {
      clab::fail("a longitude or latitude given with the points is out of "
                 "range.");
    }
  }
  if (clab_eqproj(px, py, n, p.lon0, p.lat0, earth, px, py) != CLAB_OK) {
    clab::fail("a point given with the points lies at the antipode of "
               "their centre.");
  }
}

//' @srrstats {SP6.0} back() returns centres to longitude and latitude
//'   through clab_eqinv, whose round trip with clab_eqproj is tested in
//'   src/t_geo2.cpp and test-sp.R.
//' @noRd
void back(const pts &p, double *px, double *py, int n) {
  if (!p.ll) return;
  // a projected point is within 2 R of the centre by construction, and a
  // NaN one stays NaN, so the inverse cannot fail here
  clab::kchk(clab_eqinv(px, py, n, p.lon0, p.lat0, earth, px, py), nullptr);
}

//' @srrstats {G2.7} Per-unit vectors given separately (x and y, weights,
//'   times) are matched by name when every one is named, as the package's
//'   per-group arguments are (test-sp.R).
//' @noRd
int byname(SEXP *v, const char *const *nm, int k) {
  int first = -1, cnt = 0;
  R_xlen_t n = 0;
  std::vector<std::vector<std::pair<std::string, R_xlen_t>>> key(static_cast<std::size_t>(k));
  for (int j = 0; j < k; ++j) {
    SEXP x = v[j];
    if (Rf_isNull(x)) continue;
    SEXP nms = Rf_getAttrib(x, R_NamesSymbol);
    if (!Rf_isVectorAtomic(x) || !Rf_isNull(Rf_getAttrib(x, R_DimSymbol)) ||
        Rf_isNull(nms)) {
      return 0;
    }
    if (first < 0) {
      first = j;
      n = XLENGTH(x);
    } else if (XLENGTH(x) != n) {
      return 0;
    }
    ++cnt;
    auto &kv = key[static_cast<std::size_t>(j)];
    kv.reserve(static_cast<std::size_t>(n));
    for (R_xlen_t i = 0; i < n; ++i) {
      SEXP c = STRING_ELT(nms, i);
      if (c == NA_STRING) return 0;
      kv.emplace_back(Rf_translateCharUTF8(c), i);
    }
    std::sort(kv.begin(), kv.end());
    for (std::size_t i = 1; i < kv.size(); ++i) {
      if (kv[i].first == kv[i - 1].first) return 0;  // repeated labels
    }
  }
  if (cnt < 2) return 0;
  int np = 0;
  SEXP f = v[first], fn = Rf_getAttrib(f, R_NamesSymbol);
  std::vector<R_xlen_t> o(static_cast<std::size_t>(n));
  for (int j = first + 1; j < k; ++j) {
    if (Rf_isNull(v[j])) continue;
    const auto &kv = key[static_cast<std::size_t>(j)];
    bool id = true;
    for (R_xlen_t i = 0; i < n; ++i) {
      const std::string a = Rf_translateCharUTF8(STRING_ELT(fn, i));
      auto it = std::lower_bound(kv.begin(), kv.end(), std::make_pair(a, R_xlen_t{0}));
      if (it == kv.end() || it->first != a) {
        clab::fail("`%s` must name the same units as `%s`; `%s` has no %s.",
                   nm[j], nm[first], nm[j], a.c_str());
      }
      o[static_cast<std::size_t>(i)] = it->second;
      id = id && it->second == i;
    }
    if (id) continue;
    SEXP x = v[j];
    SEXP y = PROTECT(Rf_allocVector(TYPEOF(x), n));
    ++np;
    for (R_xlen_t i = 0; i < n; ++i) {
      const R_xlen_t s = o[static_cast<std::size_t>(i)];
      switch (TYPEOF(x)) {
        case REALSXP: REAL(y)[i] = REAL_RO(x)[s]; break;
        case INTSXP: INTEGER(y)[i] = INTEGER_RO(x)[s]; break;
        case LGLSXP: LOGICAL(y)[i] = LOGICAL_RO(x)[s]; break;
        case STRSXP: SET_STRING_ELT(y, i, STRING_ELT(x, s)); break;
        default: clab::fail("`%s` must be a vector of numbers or labels.", nm[j]);
      }
    }
    Rf_copyMostAttrib(x, y);
    Rf_setAttrib(y, R_NamesSymbol, fn);
    v[j] = y;
  }
  return np;
}

SEXP gvec(SEXP x, const char *nm) {
  if (TYPEOF(x) == VECSXP && !Rf_inherits(x, "sf")) {
    const R_xlen_t n = XLENGTH(x);
    SEXP v = PROTECT(Rf_allocVector(REALSXP, n));
    for (R_xlen_t i = 0; i < n; ++i) {
      SEXP e = VECTOR_ELT(x, i);
      if (!clab::isnum(e) || XLENGTH(e) != 1) {
        clab::fail("`%s` must be numbers: a vector, a list of single numbers "
                   "or a one-row data frame.", nm);
      }
      REAL(v)[i] = clab::val1(e);
    }
    Rf_setAttrib(v, R_NamesSymbol, Rf_getAttrib(x, R_NamesSymbol));
    UNPROTECT(1);
    return v;
  }
  SEXP d = Rf_getAttrib(x, R_DimSymbol);
  if (clab::isnum(x) && TYPEOF(d) == INTSXP && XLENGTH(d) == 2 &&
      (INTEGER_RO(d)[0] == 1 || INTEGER_RO(d)[1] == 1)) {
    const bool row = INTEGER_RO(d)[0] == 1 && INTEGER_RO(d)[1] != 1;
    SEXP dn = Rf_getAttrib(x, R_DimNamesSymbol);
    SEXP v = PROTECT(clab::dcopy(PROTECT(clab::dbl(x))));
    if (!Rf_isNull(dn)) Rf_setAttrib(v, R_NamesSymbol, VECTOR_ELT(dn, row ? 1 : 0));
    UNPROTECT(2);
    return v;
  }
  return x;
}

lab codes(SEXP x, const char *nm, R_xlen_t n) {
  if (Rf_xlength(x) != n) clab::fail("`%s` must have one label per value, %ld.", nm, static_cast<long>(n));
  SEXP s = PROTECT(clab::labs(x));
  if (Rf_isNull(s) || clab::anyna(s)) {
    clab::fail("`%s` must be labels (numbers, text or a factor) with no missing "
         "values.", nm);
  }
  SEXP o = PROTECT(Rf_allocVector(INTSXP, n));
  int *ord = INTEGER(o);
  for (R_xlen_t i = 0; i < n; ++i) ord[i] = static_cast<int>(i);
  std::vector<clab::Key> k(static_cast<std::size_t>(n));
  for (R_xlen_t i = 0; i < n; ++i) k[i] = clab::key(STRING_ELT(s, i));
  if (Rf_isFactor(x)) {
    const int *f = INTEGER_RO(x);
    std::stable_sort(ord, ord + n, [f](int a, int b) { return f[a] < f[b]; });
  } else {
    std::stable_sort(ord, ord + n, [&k](int a, int b) { return clab::kcmp(k[a], k[b]) < 0; });
  }
  lab r{PROTECT(Rf_allocVector(INTSXP, n)), R_NilValue, 0};
  int g = -1;
  for (R_xlen_t q = 0; q < n; ++q) {
    const int i = ord[q];
    if (q == 0 || clab::kcmp(k[ord[q - 1]], k[i]) != 0) ++g;
    INTEGER(r.code)[i] = g;
  }
  r.g = g + 1;
  r.lev = PROTECT(Rf_allocVector(STRSXP, r.g));
  for (R_xlen_t i = 0; i < n; ++i) SET_STRING_ELT(r.lev, INTEGER(r.code)[i], STRING_ELT(s, i));
  UNPROTECT(4);
  PROTECT(r.code);
  PROTECT(r.lev);
  return r;
}

// 0 or 1 outcomes: logical or numbers
SEXP bin01(SEXP y, const char *nm, R_xlen_t n) {
  SEXP v = PROTECT(Rf_allocVector(INTSXP, n));
  if (Rf_xlength(y) != n || (TYPEOF(y) != LGLSXP && !clab::isnum(y))) {
    clab::fail("`%s` must be 0 or 1 (or TRUE or FALSE), one per case.", nm);
  }
  for (R_xlen_t i = 0; i < n; ++i) {
    const double d = TYPEOF(y) == LGLSXP
                         ? (LOGICAL_RO(y)[i] == NA_LOGICAL ? NA_REAL : LOGICAL_RO(y)[i])
                         : (TYPEOF(y) == INTSXP ? (INTEGER_RO(y)[i] == NA_INTEGER ? NA_REAL : INTEGER_RO(y)[i])
                                                : REAL_RO(y)[i]);
    if (!(d == 0 || d == 1)) clab::fail("`%s` must be 0 or 1 (or TRUE or FALSE), one per case.", nm);
    INTEGER(v)[i] = static_cast<int>(d);
  }
  UNPROTECT(1);
  return v;
}

bool flag(SEXP x, const char *nm) {
  if (TYPEOF(x) != LGLSXP || XLENGTH(x) != 1 || LOGICAL_RO(x)[0] == NA_LOGICAL) {
    clab::fail("`%s` must be TRUE or FALSE.", nm);
  }
  return LOGICAL_RO(x)[0] != 0;
}

//' @srrstats {G5.9b} Functions that simulate take `seed`; NULL draws it
//'   from R's generator, so set.seed() fixes the result (test-sp.R).
//' @noRd
unsigned int seed(SEXP s) {
  if (!Rf_isNull(s)) {
    return static_cast<unsigned int>(clab::whole(s, "seed", 0));
  }
  GetRNGstate();
  const double u = unif_rand();
  PutRNGstate();
  return static_cast<unsigned int>(u * 4294967296.0);
}

int pick(SEXP s, const char *nm, const char *const *opt, int n) {
  if (TYPEOF(s) == STRSXP && XLENGTH(s) == 1 && STRING_ELT(s, 0) != NA_STRING) {
    const char *v = CHAR(STRING_ELT(s, 0));
    for (int k = 0; k < n; ++k) {
      if (std::strcmp(v, opt[k]) == 0) return k;
    }
  }
  char all[256] = "";
  for (int k = 0; k < n; ++k) {
    std::strncat(all, k ? ", \"" : "\"", sizeof all - std::strlen(all) - 1);
    std::strncat(all, opt[k], sizeof all - std::strlen(all) - 1);
    std::strncat(all, "\"", sizeof all - std::strlen(all) - 1);
  }
  clab::fail("`%s` must be one of %s.", nm, all);
}

void win(SEXP s, double *w) {
  SEXP v = PROTECT(clab::num(s, "win", 4));
  std::memcpy(w, REAL_RO(v), 4 * sizeof(double));
  UNPROTECT(1);
  if (!(w[0] < w[1] && w[2] < w[3])) {
    clab::fail("`win` must be c(xmin, xmax, ymin, ymax) with xmin < xmax "
               "and ymin < ymax.");
  }
}

}  // namespace csp
