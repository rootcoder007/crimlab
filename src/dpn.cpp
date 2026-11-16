// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Dynamic and heterogeneous panels (crimlab.h): linear GMM on a stacked
// panel (one-step and two-step, the robust and the Windmeijer-corrected
// variances, the Sargan, Hansen and difference-in-Hansen tests and the
// Arellano-Bond tests of serial correlation), Wooldridge's test of serial
// correlation in the first-differenced regression, and Pesaran's common
// correlated effects estimators. The work arrays are made here; a failed
// allocation would end the program, so every work array is counted
// first: a size past 2^26 elements is refused with CLAB_EARG before
// anything is made. Matrices are column-major.
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <utility>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"
#include "pev.h"

namespace {

using std::size_t;
using vec = std::vector<double>;

inline size_t z(int i) { return static_cast<size_t>(i); }

// the largest work array, in elements
constexpr double cap = 67108864.0;

// The inverse of the symmetric m x m matrix a into out: by Cholesky when
// a is positive definite (1), else the Moore-Penrose inverse from the
// eigendecomposition, eigenvalues up to m eps times the largest dropped
// (2); 0 when a has no positive eigenvalue or the decomposition fails.
int sinv(const double *a, int m, double *out) {
  const size_t mm = z(m) * z(m);
  vec c(a, a + mm);
  if (pev::chol(c.data(), m)) {
    pev::cinv(c.data(), m, out);
    return 1;
  }
  vec b(a, a + mm), val(z(m)), v(mm);
  if (!gla::jacobi(b.data(), m, val.data(), v.data()) || !(val[0] > 0)) return 0;
  const double tol = val[0] * m * DBL_EPSILON;
  std::fill(out, out + mm, 0.0);
  for (int e = 0; e < m && val[z(e)] > tol; ++e) {
    const double r = 1.0 / val[z(e)];
    const double *u = v.data() + z(e) * z(m);
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) out[z(j) * z(m) + z(i)] += r * u[i] * u[j];
    }
  }
  return 2;
}

// One GMM step with weight a (m x m): p = (Z'X)' a (k x m), mi = (p Z'X)^-1
// (k x k), b = mi p Z'y. False when p Z'X is singular.
struct step {
  vec p, mi, b;
};

bool gstep(const vec &zx, const vec &zy, const double *a, int m, int k, step &s) {
  s.p.assign(z(k) * z(m), 0.0);
  for (int c = 0; c < m; ++c) {
    for (int r = 0; r < k; ++r) {
      double t = 0.0;
      for (int d = 0; d < m; ++d) t += zx[z(r) * z(m) + z(d)] * a[z(c) * z(m) + z(d)];
      s.p[z(c) * z(k) + z(r)] = t;
    }
  }
  vec mm(z(k) * z(k));
  for (int q = 0; q < k; ++q) {
    for (int r = 0; r < k; ++r) {
      double t = 0.0;
      for (int c = 0; c < m; ++c) t += s.p[z(c) * z(k) + z(r)] * zx[z(q) * z(m) + z(c)];
      mm[z(q) * z(k) + z(r)] = t;
    }
  }
  if (!pev::chol(mm.data(), k)) return false;
  s.mi.assign(z(k) * z(k), 0.0);
  pev::cinv(mm.data(), k, s.mi.data());
  vec pz(z(k), 0.0);
  for (int r = 0; r < k; ++r) {
    for (int c = 0; c < m; ++c) pz[z(r)] += s.p[z(c) * z(k) + z(r)] * zy[z(c)];
  }
  s.b.assign(z(k), 0.0);
  for (int r = 0; r < k; ++r) {
    for (int q = 0; q < k; ++q) s.b[z(r)] += s.mi[z(q) * z(k) + z(r)] * pz[z(q)];
  }
  return true;
}

// c = a b for a (r x s) and b (s x t)
vec mul(const vec &a, const vec &b, int r, int s, int t) {
  vec c(z(r) * z(t), 0.0);
  for (int j = 0; j < t; ++j) {
    for (int h = 0; h < s; ++h) {
      const double v = b[z(j) * z(s) + z(h)];
      if (v == 0.0) continue;
      for (int i = 0; i < r; ++i) c[z(j) * z(r) + z(i)] += a[z(h) * z(r) + z(i)] * v;
    }
  }
  return c;
}

// a' (r x s, so the result is s x r)
vec tr(const vec &a, int r, int s) {
  vec c(z(r) * z(s));
  for (int j = 0; j < s; ++j) {
    for (int i = 0; i < r; ++i) c[z(i) * z(s) + z(j)] = a[z(j) * z(r) + z(i)];
  }
  return c;
}

// The quadratic form u' a u of the m-vector u.
double quad(const vec &u, const double *a, int m) {
  double s = 0.0;
  for (int j = 0; j < m; ++j) {
    double t = 0.0;
    for (int i = 0; i < m; ++i) t += a[z(j) * z(m) + z(i)] * u[z(i)];
    s += t * u[z(j)];
  }
  return s;
}

// Z'v (m) for the rows r0..r1-1 of z (n rows).
void ztv(const double *zm, int n, int m, int r0, int r1, const double *v, double *out) {
  for (int c = 0; c < m; ++c) {
    const double *zc = zm + z(c) * z(n);
    double s = 0.0;
    for (int r = r0; r < r1; ++r) s += zc[r] * v[r];
    out[c] = s;
  }
}

// The data of clab_gmmpnl, held together.
struct gin {
  int n, k, m, ng, nd;
  const int *go, *eq, *dgo, *dt, *sub;
  const double *y, *x, *zm, *h, *yd, *xd, *hc, *hdd;
};

// The residuals y - x b over n rows.
vec resid(const double *y, const double *x, int n, int k, const vec &b) {
  vec e(y, y + n);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < n; ++i) e[z(i)] -= x[z(j) * z(n) + z(i)] * b[z(j)];
  }
  return e;
}

// The Arellano-Bond statistic of order l (Arellano and Bond 1991, eqs
// 8-9; Roodman 2009, sec. 3.5) from the differenced residuals d, the
// estimation residuals' moments ue (m x ng) and the step (s, v). With sig
// > 0 the error variance is sig H rather than the residuals' products
// (Roodman 2009, n. 19): sig times the blocks hdd of the differenced rows
// and hc between the estimation rows and them.
double abar(const gin &g, int l, const vec &d, const vec &ue, const step &s, const vec &v,
            double sig) {
  const int k = g.k, m = g.m;
  double num = 0.0, t1 = 0.0;
  vec q(z(k), 0.0), zv(z(m), 0.0), w, hw;
  size_t oc = 0, od = 0;
  for (int a = 0; a < g.ng; ++a) {
    const int r0 = g.dgo[a], nda = g.dgo[a + 1] - r0, ni = g.go[a + 1] - g.go[a];
    w.assign(z(nda), 0.0);
    double c = 0.0;
    for (int r = r0; r < r0 + nda; ++r) {
      int p = r - 1;
      while (p >= r0 && g.dt[p] > g.dt[r] - l) --p;
      if (p < r0 || g.dt[p] != g.dt[r] - l) continue;
      w[z(r - r0)] = d[z(p)];
      c += d[z(p)] * d[z(r)];
      for (int j = 0; j < k; ++j) q[z(j)] += d[z(p)] * g.xd[z(j) * z(g.nd) + z(r)];
    }
    num += c;
    if (sig > 0.0) {
      const double *hd = g.hdd + od, *hx = g.hc + oc;
      for (int e = 0; e < nda; ++e) {
        for (int f = 0; f < nda; ++f) t1 += sig * w[z(e)] * hd[z(f) * z(nda) + z(e)] * w[z(f)];
      }
      hw.assign(z(ni), 0.0);
      for (int f = 0; f < nda; ++f) {
        for (int e = 0; e < ni; ++e) hw[z(e)] += sig * hx[z(f) * z(ni) + z(e)] * w[z(f)];
      }
      for (int i = 0; i < m; ++i) {
        const double *zc = g.zm + z(i) * z(g.n) + z(g.go[a]);
        double t = 0.0;
        for (int e = 0; e < ni; ++e) t += zc[e] * hw[z(e)];
        zv[z(i)] += t;
      }
      od += z(nda) * z(nda);
      oc += z(ni) * z(nda);
    } else {
      t1 += c * c;
      for (int i = 0; i < m; ++i) zv[z(i)] += ue[z(a) * z(m) + z(i)] * c;
    }
  }
  // q' mi p zv and q' v q
  vec pz(z(k), 0.0);
  for (int r = 0; r < k; ++r) {
    for (int c = 0; c < m; ++c) pz[z(r)] += s.p[z(c) * z(k) + z(r)] * zv[z(c)];
  }
  double t2 = 0.0, t3 = 0.0;
  for (int r = 0; r < k; ++r) {
    for (int c = 0; c < k; ++c) {
      t2 += q[z(r)] * s.mi[z(c) * z(k) + z(r)] * pz[z(c)];
      t3 += q[z(r)] * v[z(c) * z(k) + z(r)] * q[z(c)];
    }
  }
  const double var = t1 - 2.0 * t2 + t3;
  return var > 0.0 ? num / std::sqrt(var) : NAN;
}

// The overidentification statistic without the instruments flagged in
// sub, re-estimated with the submatrix of the moment variance sv (m x m)
// as the weight's inverse (Roodman 2009, n. 24); NaN when the rest do not
// identify the coefficients. *gen is set when the inverse is general.
double subj(const vec &sv, const vec &zx, const vec &zy, const int *sub, int m, int k, int *gen) {
  std::vector<int> keep;
  for (int i = 0; i < m; ++i) {
    if (!sub[i]) keep.push_back(i);
  }
  const int ms = static_cast<int>(keep.size());
  if (ms < k || ms == m) return NAN;
  vec ss(z(ms) * z(ms)), as(z(ms) * z(ms)), zxs(z(ms) * z(k)), zys(z(ms));
  for (int j = 0; j < ms; ++j) {
    for (int i = 0; i < ms; ++i) ss[z(j) * z(ms) + z(i)] = sv[z(keep[z(j)]) * z(m) + z(keep[z(i)])];
    zys[z(j)] = zy[z(keep[z(j)])];
    for (int h = 0; h < k; ++h) zxs[z(h) * z(ms) + z(j)] = zx[z(h) * z(m) + z(keep[z(j)])];
  }
  step su;
  const int c = sinv(ss.data(), ms, as.data());
  if (!c || !gstep(zxs, zys, as.data(), ms, k, su)) return NAN;
  *gen |= c == 2;
  for (int h = 0; h < k; ++h) {
    for (int i = 0; i < ms; ++i) zys[z(i)] -= zxs[z(h) * z(ms) + z(i)] * su.b[z(h)];
  }
  return quad(zys, as.data(), ms);
}

int gmm(const gin &g, int steps, int rob, int nar, double arf, double *b, double *v, double *st,
        double *ar) {
  const int n = g.n, k = g.k, m = g.m;
  const size_t mm = z(m) * z(m);
  int gen = 0;
  // Z'HZ, block by block
  vec zhz(mm, 0.0), hz;
  size_t ho = 0;
  double htr = 0.0;
  for (int a = 0; a < g.ng; ++a) {
    const int r0 = g.go[a], ni = g.go[a + 1] - r0;
    const double *ha = g.h + ho;
    if (kutil::stop(static_cast<double>(ni) * m * (ni + m))) return CLAB_EINTR;
    hz.assign(z(ni) * z(m), 0.0);
    for (int c = 0; c < m; ++c) {
      const double *zc = g.zm + z(c) * z(n) + z(r0);
      for (int s = 0; s < ni; ++s) {
        if (zc[s] == 0.0) continue;
        for (int r = 0; r < ni; ++r) hz[z(c) * z(ni) + z(r)] += ha[z(s) * z(ni) + z(r)] * zc[s];
      }
    }
    for (int d = 0; d < m; ++d) {
      for (int c = 0; c < m; ++c) {
        const double *zc = g.zm + z(c) * z(n) + z(r0);
        double t = 0.0;
        for (int r = 0; r < ni; ++r) t += zc[r] * hz[z(d) * z(ni) + z(r)];
        zhz[z(d) * z(m) + z(c)] += t;
      }
    }
    for (int r = 0; r < ni; ++r) {
      if (g.eq[r0 + r]) htr += ha[z(r) * z(ni) + z(r)];
    }
    ho += z(ni) * z(ni);
  }
  vec a1(mm), a2(mm), s2(mm, 0.0);
  int c = sinv(zhz.data(), m, a1.data());
  if (kutil::stop(0)) return CLAB_EINTR;
  if (!c) return CLAB_EDOM;
  gen |= c == 2;
  vec zx(z(m) * z(k)), zy(z(m));
  for (int j = 0; j < k; ++j) ztv(g.zm, n, m, 0, n, g.x + z(j) * z(n), zx.data() + z(j) * z(m));
  ztv(g.zm, n, m, 0, n, g.y, zy.data());
  step s1, s2s;
  if (!gstep(zx, zy, a1.data(), m, k, s1)) return CLAB_EDOM;
  const vec e1 = resid(g.y, g.x, n, k, s1.b);
  vec u(z(m)), ue1(z(m) * z(g.ng));
  for (int a = 0; a < g.ng; ++a) {
    if (kutil::stop(static_cast<double>(m) * (m + g.go[a + 1] - g.go[a]))) return CLAB_EINTR;
    ztv(g.zm, n, m, g.go[a], g.go[a + 1], e1.data(), ue1.data() + z(a) * z(m));
    const double *w = ue1.data() + z(a) * z(m);
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) s2[z(j) * z(m) + z(i)] += w[i] * w[j];
    }
  }
  c = sinv(s2.data(), m, a2.data());
  if (kutil::stop(0)) return CLAB_EINTR;
  if (!c) return CLAB_EDOM;
  gen |= c == 2;
  if (!gstep(zx, zy, a2.data(), m, k, s2s)) return CLAB_EDOM;
  const vec e2 = resid(g.y, g.x, n, k, s2s.b);
  // the robust one-step variance (Roodman 2009, eq 15)
  const vec ps = mul(s1.p, s2, k, m, m), q1 = mul(ps, tr(s1.p, k, m), k, m, k);
  const vec v1 = mul(mul(s1.mi, q1, k, k, k), s1.mi, k, k, k);
  vec vv = s2s.mi;
  if (steps == 2 && rob) {
    // Windmeijer (2005, eq 3.4; Roodman 2009, eq 18): D's column p is
    // mi p_2 sum_i (a_ip b_i' + b_i a_ip') g, a_ip = Z_i'x_p, b_i = Z_i'e1_i,
    // g = A_2 Z'e_2
    vec ze2(zy);
    for (int j = 0; j < k; ++j) {
      for (int i = 0; i < m; ++i) ze2[z(i)] -= zx[z(j) * z(m) + z(i)] * s2s.b[z(j)];
    }
    vec gg(z(m), 0.0), uu(z(m) * z(k), 0.0), cz(z(m) * z(k));
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) gg[z(i)] += a2[z(j) * z(m) + z(i)] * ze2[z(j)];
    }
    for (int a = 0; a < g.ng; ++a) {
      if (kutil::stop(static_cast<double>(m) * k * (g.go[a + 1] - g.go[a]))) return CLAB_EINTR;
      for (int j = 0; j < k; ++j) {
        ztv(g.zm, n, m, g.go[a], g.go[a + 1], g.x + z(j) * z(n), cz.data() + z(j) * z(m));
      }
      const double *w = ue1.data() + z(a) * z(m);
      double bg = 0.0;
      for (int i = 0; i < m; ++i) bg += w[i] * gg[z(i)];
      for (int j = 0; j < k; ++j) {
        double ag = 0.0;
        for (int i = 0; i < m; ++i) ag += cz[z(j) * z(m) + z(i)] * gg[z(i)];
        for (int i = 0; i < m; ++i) {
          uu[z(j) * z(m) + z(i)] += cz[z(j) * z(m) + z(i)] * bg + w[i] * ag;
        }
      }
    }
    const vec dd = mul(mul(s2s.mi, s2s.p, k, k, m), uu, k, m, k), dt = tr(dd, k, k);
    const vec a = mul(dd, s2s.mi, k, k, k), bb = mul(s2s.mi, dt, k, k, k);
    const vec cc = mul(mul(dd, v1, k, k, k), dt, k, k, k);
    for (size_t q = 0; q < z(k) * z(k); ++q) vv[q] += a[q] + bb[q] + cc[q];
  }
  // Sargan (Arellano and Bond 1991, p. 282) and Hansen (Roodman 2009, sec. 2.5)
  vec z1(zy), z2(zy);
  for (int j = 0; j < k; ++j) {
    for (int i = 0; i < m; ++i) {
      z1[z(i)] -= zx[z(j) * z(m) + z(i)] * s1.b[z(j)];
      z2[z(i)] -= zx[z(j) * z(m) + z(i)] * s2s.b[z(j)];
    }
  }
  double se = 0.0;
  for (int r = 0; r < n; ++r) {
    if (g.eq[r]) se += e1[z(r)] * e1[z(r)];
  }
  const double sig = htr > 0.0 ? se / htr : NAN;
  // the one-step variance: robust (Roodman 2009, eq 15), or sig (X'Z A Z'X)^-1,
  // eq 14 with Var(z e) = sig Z'HZ
  if (steps == 1) {
    vv = v1;
    if (!rob) {
      for (size_t q = 0; q < z(k) * z(k); ++q) vv[q] = sig * s1.mi[q];
    }
  }
  const step &sc = steps == 1 ? s1 : s2s;
  for (int j = 0; j < k; ++j) b[j] = sc.b[z(j)];
  for (size_t q = 0; q < z(k) * z(k); ++q) v[q] = vv[q];
  st[0] = sig > 0.0 ? quad(z1, a1.data(), m) / sig : NAN;
  st[1] = quad(z2, a2.data(), m);
  st[2] = st[5] = NAN;
  st[3] = sig;
  if (g.sub) {
    st[2] = subj(s2, zx, zy, g.sub, m, k, &gen);
    if (kutil::stop(0)) return CLAB_EINTR;
    st[5] = subj(zhz, zx, zy, g.sub, m, k, &gen) / sig;
    if (kutil::stop(0)) return CLAB_EINTR;
  }
  st[4] = gen;
  // the Arellano-Bond tests on the differenced residuals of the step used
  const vec &ec = steps == 1 ? e1 : e2;
  vec ue(z(m) * z(g.ng));
  for (int a = 0; a < g.ng; ++a) {
    ztv(g.zm, n, m, g.go[a], g.go[a + 1], ec.data(), ue.data() + z(a) * z(m));
  }
  const vec d = resid(g.yd, g.xd, g.nd, k, sc.b);
  for (int l = 1; l <= nar; ++l) {
    if (kutil::stop(static_cast<double>(g.nd) * (k + l) + static_cast<double>(m) * g.ng)) {
      return CLAB_EINTR;
    }
    ar[l - 1] = abar(g, l, d, ue, sc, vv, steps == 1 && !rob ? sig * arf : 0.0);
  }
  return CLAB_OK;
}

}  // namespace

extern "C" int clab_gmmpnl(int n, int k, int m, int ng, const int *go, const double *y,
                           const double *x, const double *zm, const double *h, const int *eq,
                           int nd, const int *dgo, const int *dt, const double *yd,
                           const double *xd, const double *hc, const double *hdd, const int *sub,
                           int steps, int rob, int nar, double arf, double *b, double *v,
                           double *st, double *ar) CLAB_TRY {
  if (n < 1 || k < 1 || m < k || ng < 1 || nd < 0 || nar < 0 || (steps != 1 && steps != 2) ||
      !(arf > 0.0) || !std::isfinite(arf)) {
    return CLAB_EARG;
  }
  if (go[0] != 0 || go[ng] != n || dgo[0] != 0 || dgo[ng] != nd) return CLAB_EARG;
  for (int a = 0; a < ng; ++a) {
    if (go[a + 1] < go[a] || dgo[a + 1] < dgo[a]) return CLAB_EARG;
    for (int r = dgo[a] + 1; r < dgo[a + 1]; ++r) {
      if (dt[r] <= dt[r - 1]) return CLAB_EARG;
    }
  }
  for (size_t q = 0; q < z(n) * z(k); ++q) {
    if (!std::isfinite(x[q])) return CLAB_EDOM;
  }
  for (size_t q = 0; q < z(n) * z(m); ++q) {
    if (!std::isfinite(zm[q])) return CLAB_EDOM;
  }
  // the work arrays: m x m weights, the groups' moments (m x ng), a
  // group's H Z (n_g x m), Z'X (m x k)
  int big = 0;
  for (int a = 0; a < ng; ++a) big = std::max(big, go[a + 1] - go[a]);
  const double dm = m;
  if (dm * dm > cap || dm * ng > cap || dm * big > cap || dm * k > cap ||
      static_cast<double>(k) * k > cap) {
    return CLAB_EARG;
  }
  const gin g{n, k, m, ng, nd, go, eq, dgo, dt, sub, y, x, zm, h, yd, xd, hc, hdd};
  return gmm(g, steps, rob, nar, arf, b, v, st, ar);
} CLAB_CATCH(CLAB_ENOMEM)

namespace {

int wsct(int n, int p, const double *y, const double *x, const int *u, const int *t, int cons,
         int raw, double *b, double *vb, double *out) {
  std::vector<int> o(z(n));
  std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](int a, int c) {
    return u[a] != u[c] ? u[a] < u[c] : t[a] < t[c];
  });
  // the first differences: rows whose unit has the period before
  std::vector<int> fr, fu, ft;
  for (int q = 1; q < n; ++q) {
    const int r = o[z(q)], s = o[z(q - 1)];
    if (u[r] == u[s] && t[r] == t[s] + 1) {
      fr.push_back(q);
      fu.push_back(u[r]);
      ft.push_back(t[r]);
    }
  }
  const int nf = static_cast<int>(fr.size());
  if (nf <= p) return CLAB_EDOM;
  vec dy(z(nf)), dx(z(nf) * z(p)), e(z(nf));
  for (int i = 0; i < nf; ++i) {
    const int r = o[z(fr[z(i)])], s = o[z(fr[z(i)] - 1)];
    dy[z(i)] = y[r] - y[s];
    for (int j = 0; j < p; ++j) {
      dx[z(j) * z(nf) + z(i)] = x[z(j) * z(n) + z(r)] - x[z(j) * z(n) + z(s)];
    }
  }
  // clusters numbered from 0 among the rows used
  auto recode = [](const std::vector<int> &c, std::vector<int> &out) {
    std::vector<int> l(c);
    std::sort(l.begin(), l.end());
    l.erase(std::unique(l.begin(), l.end()), l.end());
    out.resize(c.size());
    for (size_t i = 0; i < c.size(); ++i) {
      out[i] = static_cast<int>(std::lower_bound(l.begin(), l.end(), c[i]) - l.begin());
    }
    return static_cast<int>(l.size());
  };
  auto fit = [&](const vec &xx, const vec &yy, int nn, int kk, const std::vector<int> &cl, int g,
                 double *bb, double *vv, vec &r) {
    vec xi(z(kk) * z(kk)), wk(z(kk) * z(kk) + z(kk)), wc(z(g) * z(kk) + 2 * z(kk) * z(kk));
    r.assign(z(nn), 0.0);
    if (pev::ols(xx.data(), yy.data(), nn, kk, nullptr, bb, xi.data(), r.data(), wk.data()) !=
        CLAB_OK) {
      return false;
    }
    pev::crv(xx.data(), r.data(), nullptr, nn, kk, xi.data(), cl.data(), g, nullptr, 0, nullptr, 0,
             kk, vv, wc.data());
    if (raw) {
      const double f = static_cast<double>(g) / (g - 1) * (nn - 1.0) / (nn - kk);
      for (int q = 0; q < kk * kk; ++q) vv[q] /= f;
    }
    return true;
  };
  std::vector<int> cl;
  const int g1 = recode(fu, cl);
  if (p > 0) {
    if (g1 < 2 || !fit(dx, dy, nf, p, cl, g1, b, vb, e)) return CLAB_EDOM;
  } else {
    e = dy;
  }
  // the residual on its lag: pairs of differences in consecutive periods
  std::vector<int> pu;
  vec pe, pl;
  for (int i = 1; i < nf; ++i) {
    if (fu[z(i)] == fu[z(i - 1)] && ft[z(i)] == ft[z(i - 1)] + 1) {
      pu.push_back(fu[z(i)]);
      pe.push_back(e[z(i)]);
      pl.push_back(e[z(i - 1)]);
    }
  }
  const int np = static_cast<int>(pe.size()), k2 = cons ? 2 : 1;
  std::vector<int> c2;
  const int g2 = recode(pu, c2);
  if (np <= k2 || g2 < 2) return CLAB_EDOM;
  vec x2(z(np) * z(k2), 1.0), r2;
  std::copy(pl.begin(), pl.end(), x2.begin());
  double b2[2], v2[4];
  if (!fit(x2, pe, np, k2, c2, g2, b2, v2, r2) || !(v2[0] > 0.0)) return CLAB_EDOM;
  const double se = std::sqrt(v2[0]), f = (b2[0] + 0.5) / se;
  out[0] = b2[0];
  out[1] = se;
  out[2] = f * f;
  out[3] = g2;
  out[4] = np;
  out[5] = nf;
  return CLAB_OK;
}

// What the dynamic estimator adds to Pesaran's: p lags of the averages
// kept, the rows of each unit's regression (rm: nu x nt flags, unit i at
// rm[i nt + t]; NULL for every observed row), whether the pooled
// estimator is made, and the heterogeneity-bias statistic (hb: 2).
struct cex {
  int p = 0;
  const int *rm = nullptr;
  int pool = 1;
  double *hb = nullptr;
};

// The heterogeneity-bias statistic of Juodis and Reese (2026, case study
// 2): delta_i = (A_i - mean_j A_j) (b_i - b_MG), A_i = X_i'M_i X_i, and
// HB = delta' V^-1 delta with delta the mean of the delta_i and V = sum
// delta_i delta_i' / N^2, chi-squared on k degrees of freedom; NaN when V
// is singular.
double hbias(int nu, int k, const vec &ai, const double *bi, const double *bmg) {
  const size_t kk = z(k) * z(k);
  vec am(kk, 0.0), dl(z(nu) * z(k)), dm(z(k), 0.0), vv(kk, 0.0), vi(kk);
  for (int i = 0; i < nu; ++i) {
    for (size_t q2 = 0; q2 < kk; ++q2) am[q2] += ai[z(i) * kk + q2] / nu;
  }
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int g = 0; g < k; ++g) {
        s += (ai[z(i) * kk + z(g) * z(k) + z(j)] - am[z(g) * z(k) + z(j)]) *
             (bi[z(i) * z(k) + z(g)] - bmg[g]);
      }
      dl[z(i) * z(k) + z(j)] = s;
      dm[z(j)] += s / nu;
    }
  }
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) {
        vv[z(g) * z(k) + z(j)] +=
            dl[z(i) * z(k) + z(j)] * dl[z(i) * z(k) + z(g)] / (static_cast<double>(nu) * nu);
      }
    }
  }
  if (!pev::chol(vv.data(), k)) return NAN;
  pev::cinv(vv.data(), k, vi.data());
  return quad(dm, vi.data(), k);
}

// Pesaran (2006): the unit estimates, CCEMG and CCEP with their
// variances. A NaN in y marks a unit-period not observed: the averages of
// a period are over the units observed in it, their weights scaled to sum
// to 1 there, and each unit's regression is over its own periods. With
// ex.p > 0 the averages kept enter H with their lags 1..p too (Chudik and
// Pesaran 2015, eq 24), and with ex.rm each unit's regression is over the
// rows flagged; only the columns of x whose averages are kept enter the
// averages.
int cce(int nu, int nt, int k, const double *y, const double *x, const double *w, const double *d,
        int q, const int *av, const cex &ex, double *bi, double *bmg, double *vmg, double *bp,
        double *vp, int *bad) {
  // H = (D, ybar_w, xbar_w), T x (q + 1 + k), the averages flagged in av kept
  const int ha = q + 1 + k;
  vec hm(z(nt) * z(ha), 0.0), sw(z(nt), 0.0);
  std::vector<int> cnt(z(nt), 0);
  for (int c = 0; c < q; ++c) {
    for (int t = 0; t < nt; ++t) hm[z(c) * z(nt) + z(t)] = d[z(c) * z(nt) + z(t)];
  }
  for (int i = 0; i < nu; ++i) {
    for (int t = 0; t < nt; ++t) {
      if (std::isnan(y[z(i) * z(nt) + z(t)])) continue;
      sw[z(t)] += w[i];
      ++cnt[z(t)];
      hm[z(q) * z(nt) + z(t)] += w[i] * y[z(i) * z(nt) + z(t)];
      for (int j = 0; j < k; ++j) {
        if (!av[j + 1]) continue;
        hm[z(q + 1 + j) * z(nt) + z(t)] += w[i] * x[(z(j) * z(nu) + z(i)) * z(nt) + z(t)];
      }
    }
  }
  for (int t = 0; t < nt; ++t) {
    if (cnt[z(t)] < 2) {
      *bad = -2 - t;
      return CLAB_EDOM;
    }
    // a period every unit has keeps the weights as they are
    if (cnt[z(t)] == nu) continue;
    for (int c = q; c < ha; ++c) hm[z(c) * z(nt) + z(t)] /= sw[z(t)];
  }
  int hc = q;
  for (int c = 0; c <= k; ++c) {
    if (!av[c]) continue;
    if (q + c != hc) {
      std::copy(hm.begin() + static_cast<std::ptrdiff_t>(z(q + c) * z(nt)),
                hm.begin() + static_cast<std::ptrdiff_t>(z(q + c + 1) * z(nt)),
                hm.begin() + static_cast<std::ptrdiff_t>(z(hc) * z(nt)));
    }
    ++hc;
  }
  if (hc == 0) return CLAB_EARG;
  // the lags 1..p of the averages kept, after them (NaN before period 0;
  // the rows flagged start at period p)
  if (ex.p > 0) {
    const int na = hc - q, hl = hc + na * ex.p;
    if (z(nt) * z(hl) > hm.size()) hm.resize(z(nt) * z(hl), NAN);
    for (int l = 1; l <= ex.p; ++l) {
      for (int c = 0; c < na; ++c) {
        double *o = hm.data() + z(hc + (l - 1) * na + c) * z(nt);
        const double *s = hm.data() + z(q + c) * z(nt);
        for (int t = 0; t < l; ++t) o[t] = NAN;
        for (int t = l; t < nt; ++t) o[t] = s[t - l];
      }
    }
    hc = hl;
  }
  const size_t kk = z(k) * z(k);
  vec ai(z(nu) * kk), ci(z(nu) * z(k)), ti(z(nu)), l(kk), li(kk);
  std::vector<int> per;
  for (int i = 0; i < nu; ++i) {
    per.clear();
    for (int t = 0; t < nt; ++t) {
      const bool on = ex.rm ? ex.rm[z(i) * z(nt) + z(t)] != 0 : !std::isnan(y[z(i) * z(nt) + z(t)]);
      if (on) per.push_back(t);
    }
    const int ni = static_cast<int>(per.size());
    ti[z(i)] = ni;
    const double wk =
        static_cast<double>(ni) * (k + 1) * (hc + k) + static_cast<double>(hc) * hc * hc;
    if (kutil::stop(wk)) {
      return CLAB_EINTR;
    }
    // the unit's rows of H, (H_i'H_i)^-, and M_i v = v - H_i (H_i'H_i)^- H_i'v
    vec h(z(ni) * z(hc)), hh(z(hc) * z(hc)), hi(z(hc) * z(hc)), hv(z(hc)), gv(z(hc));
    for (int a = 0; a < hc; ++a) {
      for (int r = 0; r < ni; ++r) h[z(a) * z(ni) + z(r)] = hm[z(a) * z(nt) + z(per[z(r)])];
    }
    for (int a = 0; a < hc; ++a) {
      for (int c = 0; c < hc; ++c) {
        double s = 0.0;
        for (int r = 0; r < ni; ++r) s += h[z(a) * z(ni) + z(r)] * h[z(c) * z(ni) + z(r)];
        hh[z(c) * z(hc) + z(a)] = s;
      }
    }
    if (!sinv(hh.data(), hc, hi.data())) {
      *bad = i;
      return CLAB_EDOM;
    }
    auto proj = [&](double *v) {
      for (int a = 0; a < hc; ++a) {
        double s = 0.0;
        for (int r = 0; r < ni; ++r) s += h[z(a) * z(ni) + z(r)] * v[r];
        hv[z(a)] = s;
      }
      for (int a = 0; a < hc; ++a) {
        double s = 0.0;
        for (int c = 0; c < hc; ++c) s += hi[z(c) * z(hc) + z(a)] * hv[z(c)];
        gv[z(a)] = s;
      }
      for (int r = 0; r < ni; ++r) {
        double s = 0.0;
        for (int a = 0; a < hc; ++a) s += h[z(a) * z(ni) + z(r)] * gv[z(a)];
        v[r] -= s;
      }
    };
    vec mx(z(ni) * z(k)), my(z(ni)), rn(z(k), 0.0);
    for (int j = 0; j < k; ++j) {
      for (int r = 0; r < ni; ++r) {
        const double e = x[(z(j) * z(nu) + z(i)) * z(nt) + z(per[z(r)])];
        mx[z(j) * z(ni) + z(r)] = e;
        rn[z(j)] += e * e;
      }
      proj(mx.data() + z(j) * z(ni));
    }
    for (int r = 0; r < ni; ++r) my[z(r)] = y[z(i) * z(nt) + z(per[z(r)])];
    proj(my.data());
    double *a = ai.data() + z(i) * kk, *c = ci.data() + z(i) * z(k);
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int r = 0; r < ni; ++r) s += mx[z(j) * z(ni) + z(r)] * my[z(r)];
      c[j] = s;
      for (int g = 0; g < k; ++g) {
        double e = 0.0;
        for (int r = 0; r < ni; ++r) e += mx[z(j) * z(ni) + z(r)] * mx[z(g) * z(ni) + z(r)];
        a[z(g) * z(k) + z(j)] = e;
      }
    }
    // a regressor the averages take out (to 1e-7 of its norm) leaves only
    // rounding: Cholesky's relative pivot cannot see that in one column
    bool gone = false;
    for (int j = 0; j < k; ++j) gone = gone || !(a[z(j) * z(k) + z(j)] > 1e-14 * rn[z(j)]);
    std::copy(a, a + kk, l.begin());
    if (gone || !pev::chol(l.data(), k)) {
      *bad = i;
      return CLAB_EDOM;
    }
    pev::cinv(l.data(), k, li.data());
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int g = 0; g < k; ++g) s += li[z(g) * z(k) + z(j)] * c[g];
      bi[z(i) * z(k) + z(j)] = s;
    }
  }
  // CCEMG (eqs 53 and 58)
  std::fill(bmg, bmg + k, 0.0);
  std::fill(vmg, vmg + kk, 0.0);
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) bmg[j] += bi[z(i) * z(k) + z(j)] / nu;
  }
  for (int i = 0; i < nu; ++i) {
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) {
        vmg[z(g) * z(k) + z(j)] += (bi[z(i) * z(k) + z(j)] - bmg[j]) *
                                    (bi[z(i) * z(k) + z(g)] - bmg[g]) / (nu * (nu - 1.0));
      }
    }
  }
  if (ex.hb) {
    ex.hb[0] = NAN;
    ex.hb[1] = k;
  }
  if (!ex.pool) return CLAB_OK;
  if (ex.hb) ex.hb[0] = hbias(nu, k, ai, bi, bmg);
  // CCEP (eq 65, theta_i = w_i) and its variance (eqs 67-69), each unit's
  // X_i'M X_i over its own number of periods
  vec pa(kk, 0.0), pc(z(k), 0.0), ps(kk, 0.0), rs(kk, 0.0), dv(z(k)), ad(z(k));
  double w2 = 0.0;
  for (int i = 0; i < nu; ++i) {
    w2 += w[i] * w[i];
    for (size_t q2 = 0; q2 < kk; ++q2) {
      pa[q2] += w[i] * ai[z(i) * kk + q2];
      ps[q2] += w[i] * ai[z(i) * kk + q2] / ti[z(i)];
    }
    for (int j = 0; j < k; ++j) pc[z(j)] += w[i] * ci[z(i) * z(k) + z(j)];
  }
  for (int i = 0; i < nu; ++i) {
    const double wt2 = w[i] * w[i] / (w2 / nu);
    for (int j = 0; j < k; ++j) dv[z(j)] = bi[z(i) * z(k) + z(j)] - bmg[j];
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int g = 0; g < k; ++g) s += ai[z(i) * kk + z(g) * z(k) + z(j)] / ti[z(i)] * dv[z(g)];
      ad[z(j)] = s;
    }
    for (int j = 0; j < k; ++j) {
      for (int g = 0; g < k; ++g) rs[z(g) * z(k) + z(j)] += wt2 * ad[z(j)] * ad[z(g)] / (nu - 1.0);
    }
  }
  vec pl(pa), pi(kk);
  if (!pev::chol(pl.data(), k)) {
    *bad = -1;
    return CLAB_EDOM;
  }
  pev::cinv(pl.data(), k, pi.data());
  for (int j = 0; j < k; ++j) {
    double s = 0.0;
    for (int g = 0; g < k; ++g) s += pi[z(g) * z(k) + z(j)] * pc[z(g)];
    bp[j] = s;
  }
  vec ql(ps), qi(kk);
  if (!pev::chol(ql.data(), k)) {
    *bad = -1;
    return CLAB_EDOM;
  }
  pev::cinv(ql.data(), k, qi.data());
  const vec vr = mul(mul(qi, rs, k, k, k), qi, k, k, k);
  for (size_t q2 = 0; q2 < kk; ++q2) vp[q2] = w2 * vr[q2];
  return CLAB_OK;
}

}  // namespace

extern "C" int clab_wsctst(int n, int p, const double *y, const double *x, const int *u,
                           const int *t, int cons, int raw, double *b, double *vb,
                           double *out) CLAB_TRY {
  if (n < 2 || p < 0) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) return CLAB_EDOM;
  }
  for (size_t q = 0; q < z(n) * z(p); ++q) {
    if (!std::isfinite(x[q])) return CLAB_EDOM;
  }
  // the work arrays: the differences (n x p) and p x p matrices
  if (static_cast<double>(p) * p > cap || static_cast<double>(n) * (p + 4) > cap) return CLAB_EARG;
  return wsct(n, p, y, x, u, t, cons, raw, b, vb, out);
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_ccepnl(int nu, int nt, int k, const double *y, const double *x,
                           const double *w, const double *d, int q, const int *av, double *bi,
                           double *bmg, double *vmg, double *bp, double *vp, int *bad) CLAB_TRY {
  *bad = -1;
  if (nu < 2 || nt < 2 || k < 1 || q < 0) return CLAB_EARG;
  double sw = 0.0;
  for (int i = 0; i < nu; ++i) {
    if (!(w[i] > 0.0) || !std::isfinite(w[i])) return CLAB_EDOM;
    sw += w[i];
  }
  if (std::fabs(sw - 1.0) > 1e-12) return CLAB_EDOM;
  const size_t ny = z(nu) * z(nt);
  for (size_t r = 0; r < ny; ++r) {
    if (std::isnan(y[r])) continue;
    if (!std::isfinite(y[r])) return CLAB_EDOM;
    for (int j = 0; j < k; ++j) {
      if (!std::isfinite(x[z(j) * ny + r])) return CLAB_EDOM;
    }
  }
  for (size_t r = 0; r < z(nt) * z(q); ++r) {
    if (!std::isfinite(d[r])) return CLAB_EDOM;
  }
  // the work arrays: the units' k x k matrices and the averages
  if (static_cast<double>(nu) * k * k > cap || static_cast<double>(nt) * (q + k + 1) > cap) {
    return CLAB_EARG;
  }
  return cce(nu, nt, k, y, x, w, d, q, av, cex(), bi, bmg, vmg, bp, vp, bad);
} CLAB_CATCH(CLAB_ENOMEM)

extern "C" int clab_ccedyn(int nu, int nt, int k, const double *y, const double *x,
                           const double *w, const double *d, int q, const int *av, int p,
                           const int *rm, int pool, double *bi, double *bmg, double *vmg,
                           double *bp, double *vp, double *hb, int *bad) CLAB_TRY {
  *bad = -1;
  if (nu < 2 || nt < 2 || k < 1 || q < 0 || p < 0 || p >= nt || !rm) return CLAB_EARG;
  double sw = 0.0;
  for (int i = 0; i < nu; ++i) {
    if (!(w[i] > 0.0) || !std::isfinite(w[i])) return CLAB_EDOM;
    sw += w[i];
  }
  if (std::fabs(sw - 1.0) > 1e-12) return CLAB_EDOM;
  const size_t ny = z(nu) * z(nt);
  for (size_t r = 0; r < ny; ++r) {
    const bool in = rm[r] != 0;
    if (in && (std::isnan(y[r]) || static_cast<int>(r % z(nt)) < p)) return CLAB_EARG;
    if (std::isnan(y[r])) continue;
    if (!std::isfinite(y[r])) return CLAB_EDOM;
    for (int j = 0; j < k; ++j) {
      if ((in || av[j + 1]) && !std::isfinite(x[z(j) * ny + r])) return CLAB_EDOM;
    }
  }
  for (size_t r = 0; r < z(nt) * z(q); ++r) {
    if (!std::isfinite(d[r])) return CLAB_EDOM;
  }
  // the work arrays: the units' k x k matrices and the averages with lags
  const double hl = static_cast<double>(q) + (k + 1.0) * (p + 1.0);
  const double cl = kutil::maxcell;
  if (static_cast<double>(nu) * k * k > cl || static_cast<double>(nt) * hl > cl || hl * hl > cl) {
    return CLAB_EARG;
  }
  cex ex;
  ex.p = p;
  ex.rm = rm;
  ex.pool = pool != 0;
  ex.hb = hb;
  if (!ex.pool) {
    std::fill(bp, bp + k, NAN);
    std::fill(vp, vp + z(k) * z(k), NAN);
  }
  return cce(nu, nt, k, y, x, w, d, q, av, ex, bi, bmg, vmg, bp, vp, bad);
} CLAB_CATCH(CLAB_ENOMEM)
