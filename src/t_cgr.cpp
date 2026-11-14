// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Standalone tests of the centrographic kernels (src/cgr.cpp).
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include "cgeo.h"
#include "crimlab.h"
#include "tchk.h"

namespace {

// a message of the extended run (standalone, "ext"); silent by default
template <class... A>
void info(int big, const char *f, A... a) {
  if (big) std::printf(f, a...);
}

const double PI = cgeo::pi;
double gr(int i, int a, int m) {
  return static_cast<double>((i * a) % m) / m;
}
double gauss(cgeo::rng &r) {
  double u = r.unif();
  while (u == 0.0) u = r.unif();
  return std::sqrt(-2.0 * std::log(u)) * std::cos(2.0 * PI * r.unif());
}

// set A of dev/ref/cpp/cgr.R
const int NA = 200;
double ax[NA], ay[NA], aw[NA];
void setA() {
  for (int i = 0; i < NA; ++i) {
    ax[i] = gr(i, 7919, 1009) * 1000;
    ay[i] = gr(i, 104729, 1013) * 1000;
    aw[i] = 1 + (i % 3);
  }
}

// ---- exact orientation by integers, for points on the grid 2^-53 Z^2
// with |coordinate| < 2^5: differences fit in 59 bits, products in 118.
struct u128 {
  std::uint64_t hi, lo;
};
u128 mulu(std::uint64_t a, std::uint64_t b) {
  const std::uint64_t a0 = a & 0xffffffffULL, a1 = a >> 32;
  const std::uint64_t b0 = b & 0xffffffffULL, b1 = b >> 32;
  const std::uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0,
                      p11 = a1 * b1;
  const std::uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffULL) +
                            (p10 & 0xffffffffULL);
  u128 r;
  r.lo = (p00 & 0xffffffffULL) | (mid << 32);
  r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
  return r;
}
int cmpu(u128 a, u128 b) {
  if (a.hi != b.hi) return a.hi < b.hi ? -1 : 1;
  if (a.lo != b.lo) return a.lo < b.lo ? -1 : 1;
  return 0;
}
// sign of p * q - r * s for signed 64-bit p, q, r, s
int detsign(std::int64_t p, std::int64_t q, std::int64_t r, std::int64_t s) {
  const int sl = (p == 0 || q == 0) ? 0 : ((p < 0) != (q < 0) ? -1 : 1);
  const int sr = (r == 0 || s == 0) ? 0 : ((r < 0) != (s < 0) ? -1 : 1);
  auto ab = [](std::int64_t v) {
    return static_cast<std::uint64_t>(v < 0 ? -v : v);
  };
  if (sl != sr) return sl > sr ? 1 : -1;  // sign(L - R) when signs differ
  if (sl == 0) return 0;
  const int c = cmpu(mulu(ab(p), ab(q)), mulu(ab(r), ab(s)));
  return sl > 0 ? c : -c;
}
std::int64_t fix(double v) { return static_cast<std::int64_t>(v * 0x1p53); }
int iorient(double ax_, double ay_, double bx, double by, double cx,
            double cy) {
  const std::int64_t a0 = fix(ax_), a1 = fix(ay_), b0 = fix(bx), b1 = fix(by),
                     c0 = fix(cx), c1 = fix(cy);
  return detsign(a0 - c0, b1 - c1, a1 - c1, b0 - c0);
}
int norient(double ax_, double ay_, double bx, double by, double cx,
            double cy) {
  const double d = (ax_ - cx) * (by - cy) - (ay_ - cy) * (bx - cx);
  return d > 0 ? 1 : (d < 0 ? -1 : 0);
}

// brute-force hull of distinct points by gift wrapping under predicate o:
// vertex positions anticlockwise from the lowest-x (then y) point
int giftwrap(const double *x, const double *y, int n,
             int (*o)(double, double, double, double, double, double),
             int *out) {
  int s = 0;
  for (int i = 1; i < n; ++i)
    if (x[i] < x[s] || (x[i] == x[s] && y[i] < y[s])) s = i;
  int h = 0, cur = s;
  do {
    out[h++] = cur;
    int c = (cur + 1) % n;
    for (int p = 0; p < n; ++p) {
      if (p == cur) continue;
      const int t = o(x[cur], y[cur], x[c], y[c], x[p], y[p]);
      const double dc = std::hypot(x[c] - x[cur], y[c] - y[cur]);
      const double dp = std::hypot(x[p] - x[cur], y[p] - y[cur]);
      if (t < 0 || (t == 0 && dp > dc)) c = p;
    }
    cur = c;
  } while (cur != s && h < n);
  return h;
}

// compass search for the Weber point: an independent minimiser
double wobj(const double *x, const double *y, const double *w, int n,
            double cx, double cy) {
  double f = 0;
  for (int i = 0; i < n; ++i) f += (w ? w[i] : 1.0) * std::hypot(x[i] - cx, y[i] - cy);
  return f;
}
void compass(const double *x, const double *y, const double *w, int n,
             double *cx, double *cy, double step) {
  double f = wobj(x, y, w, n, *cx, *cy);
  while (step > 1e-13 * (1 + std::fabs(*cx) + std::fabs(*cy))) {
    bool moved = false;
    for (int d = 0; d < 8; ++d) {
      const double a = d * PI / 4;
      const double tx = *cx + step * std::cos(a), ty = *cy + step * std::sin(a);
      const double ft = wobj(x, y, w, n, tx, ty);
      if (ft < f) {
        f = ft;
        *cx = tx;
        *cy = ty;
        moved = true;
        break;
      }
    }
    if (!moved) step *= 0.5;
  }
}

void t_cgstat(tk &T, int /*big*/) {
  // 3 x 3 lattice, spacing 2, at (10, 20)
  double x[9], y[9], o[10];
  for (int i = 0; i < 9; ++i) {
    x[i] = 10 + 2 * (i % 3);
    y[i] = 20 + 2 * (i / 3);
  }
  tk_true(T, "cgs.ok", clab_cgstat(x, y, nullptr, 9, 1, o) == CLAB_OK);
  tk_eq(T, "cgs.W", o[0], 9, 0);
  tk_eq(T, "cgs.mx", o[1], 12, 1e-15);
  tk_eq(T, "cgs.my", o[2], 22, 1e-15);
  tk_eq(T, "cgs.sx1", o[3], std::sqrt(3.0), 1e-15);
  tk_eq(T, "cgs.sdd1", o[5], std::sqrt(6.0), 1e-15);
  clab_cgstat(x, y, nullptr, 9, 0, o);
  tk_eq(T, "cgs.sx0", o[4], std::sqrt(24.0 / 9), 1e-15);
  clab_cgstat(x, y, nullptr, 9, 2, o);
  tk_eq(T, "cgs.sdd2", o[5], std::sqrt(48.0 / 7), 1e-15);
  // 12 points on a circle of radius 5: sdd (k = 0) is the radius
  double cx[12], cy[12];
  for (int i = 0; i < 12; ++i) {
    cx[i] = 100 + 5 * std::cos(2 * PI * i / 12);
    cy[i] = -50 + 5 * std::sin(2 * PI * i / 12);
  }
  clab_cgstat(cx, cy, nullptr, 12, 0, o);
  tk_eq(T, "cgs.circle", o[5], 5, 1e-14);
  tk_eq(T, "cgs.circx", o[3], 5 / std::sqrt(2.0), 1e-14);
  tk_eq(T, "cgs.circmx", o[1], 100, 1e-15);
  // geometric and harmonic means
  double gx[3] = {1, 2, 4}, gy[3] = {3, 0, 3};
  clab_cgstat(gx, gy, nullptr, 3, 0, o);
  tk_eq(T, "cgs.gm", o[6], 2, 1e-15);
  tk_eq(T, "cgs.hm", o[8], 12.0 / 7, 1e-15);
  tk_true(T, "cgs.gmnan", std::isnan(o[7]) && std::isnan(o[9]));
  double gw[3] = {1, 0, 1};  // the zero has no weight: defined again
  clab_cgstat(gx, gy, gw, 3, 0, o);
  tk_eq(T, "cgs.gm0w", o[7], 3, 1e-15);
  tk_eq(T, "cgs.gmw", o[6], 2, 1e-15);
  // integer weights = replicated points
  setA();
  static double rx[3 * NA], ry[3 * NA];
  int m = 0;
  for (int i = 0; i < NA; ++i)
    for (int k = 0; k < aw[i]; ++k) {
      rx[m] = ax[i];
      ry[m++] = ay[i];
    }
  for (int k = 0; k <= 2; ++k) {
    double a[10], b[10];
    clab_cgstat(ax, ay, aw, NA, k, a);
    clab_cgstat(rx, ry, nullptr, m, k, b);
    for (int j = 0; j < 10; ++j) tk_eq(T, "cgs.repl", a[j], b[j], 1e-13);
  }
  // naive two-pass formulas on random data
  cgeo::rng r(11);
  for (int rep = 0; rep < 20; ++rep) {
    double px[50], py[50];
    for (int i = 0; i < 50; ++i) {
      px[i] = r.unif() * 100;
      py[i] = r.unif() * 50 + 10;
    }
    double s1 = 0, s2 = 0;
    for (int i = 0; i < 50; ++i) {
      s1 += px[i];
      s2 += py[i];
    }
    s1 /= 50;
    s2 /= 50;
    double q1 = 0, q2 = 0;
    for (int i = 0; i < 50; ++i) {
      q1 += (px[i] - s1) * (px[i] - s1);
      q2 += (py[i] - s2) * (py[i] - s2);
    }
    clab_cgstat(px, py, nullptr, 50, 1, o);
    tk_eq(T, "cgs.naivex", o[3], std::sqrt(q1 / 49), 1e-13);
    tk_eq(T, "cgs.naived", o[5], std::sqrt((q1 + q2) / 49), 1e-13);
    // rotation by 0.7 keeps sdd; a shift of 1e7 keeps every dispersion
    double qx[50], qy[50], sx_[50], sy_[50], o2[10], o3[10];
    for (int i = 0; i < 50; ++i) {
      qx[i] = std::cos(0.7) * px[i] - std::sin(0.7) * py[i];
      qy[i] = std::sin(0.7) * px[i] + std::cos(0.7) * py[i];
      sx_[i] = px[i] + 1e7;
      sy_[i] = py[i] - 1e7;
    }
    clab_cgstat(qx, qy, nullptr, 50, 1, o2);
    tk_eq(T, "cgs.rot", o2[5], o[5], 1e-13);
    clab_cgstat(sx_, sy_, nullptr, 50, 1, o3);
    tk_eq(T, "cgs.shift", o3[3], o[3], 1e-9);
    tk_eq(T, "cgs.shiftmx", o3[1] - 1e7, o[1], 1e-8);
  }
  // guards
  double one = 1.0, nanv = std::nan(""), neg = -1;
  tk_true(T, "cgs.n0", clab_cgstat(x, y, nullptr, 0, 0, o) == CLAB_EARG);
  tk_true(T, "cgs.k", clab_cgstat(x, y, nullptr, 9, 3, o) == CLAB_EARG);
  tk_true(T, "cgs.nan", clab_cgstat(&nanv, &one, nullptr, 1, 0, o) == CLAB_EDOM);
  tk_true(T, "cgs.negw", clab_cgstat(&one, &one, &neg, 1, 0, o) == CLAB_EDOM);
  double zw[2] = {0, 0};
  tk_true(T, "cgs.zerow", clab_cgstat(x, y, zw, 2, 0, o) == CLAB_EDOM);
  tk_true(T, "cgs.dof", clab_cgstat(&one, &one, nullptr, 1, 1, o) == CLAB_EDOM);
  tk_true(T, "cgs.n1", clab_cgstat(&one, &one, nullptr, 1, 0, o) == CLAB_OK &&
                   o[3] == 0 && o[5] == 0 && o[1] == 1);
}

void t_cgmedn(tk &T, int /*big*/) {
  double x[3] = {3, 1, 2}, y[3] = {-1, 5, 0}, w[4], o[2];
  tk_true(T, "cgm.ok", clab_cgmedn(x, y, 3, w, o) == CLAB_OK);
  tk_true(T, "cgm.odd", o[0] == 2 && o[1] == 0);
  double x4[4] = {4, 1, 3, 2}, y4[4] = {0, 0, 1, 1};
  clab_cgmedn(x4, y4, 4, w, o);
  tk_true(T, "cgm.even", o[0] == 2.5 && o[1] == 0.5);
  double nanv = std::nan("");
  tk_true(T, "cgm.nan", clab_cgmedn(&nanv, y, 1, w, o) == CLAB_EDOM);
  tk_true(T, "cgm.n0", clab_cgmedn(x, y, 0, w, o) == CLAB_EARG);
}

void t_webpt(tk &T, int /*big*/) {
  double o[6];
  // square: the centre, found at the start (an exact minimiser)
  double sx[4] = {0, 2, 2, 0}, sy[4] = {0, 0, 2, 2};
  tk_true(T, "web.ok", clab_webpt(sx, sy, nullptr, 4, nullptr, 1e-12, 100, o) ==
                   CLAB_OK);
  tk_eq(T, "web.sqx", o[0], 1, 1e-15);
  tk_eq(T, "web.sqy", o[1], 1, 1e-15);
  tk_true(T, "web.sqcode", o[5] == 2);
  // a triangle with an angle above 120 degrees at (0, 0): the minimiser is
  // that vertex. Plain Weiszfeld is undefined there; the modified step is
  // not, and the vertex is certified by |R| <= eta.
  double tx[3] = {0, 1, -1}, ty[3] = {0, 0, 0.2};
  clab_webpt(tx, ty, nullptr, 3, nullptr, 1e-12, 100000, o);
  tk_true(T, "web.vertex", o[0] == 0 && o[1] == 0 && o[5] == 2);
  double st[2] = {0, 0};
  clab_webpt(tx, ty, nullptr, 3, st, 1e-15, 100, o);
  tk_true(T, "web.vertexexact", o[0] == 0 && o[1] == 0 && o[5] == 2 && o[3] == 0);
  // equilateral triangle: the centroid
  double ex[3] = {0, 1, 0.5}, ey[3] = {0, 0, std::sqrt(3.0) / 2};
  clab_webpt(ex, ey, nullptr, 3, nullptr, 1e-14, 1000, o);
  tk_eq(T, "web.eqx", o[0], 0.5, 1e-13);
  tk_eq(T, "web.eqy", o[1], std::sqrt(3.0) / 6, 1e-13);
  // collinear, odd count: the middle point
  double lx[5] = {0, 1, 5, 7, 20}, ly[5];
  for (int i = 0; i < 5; ++i) ly[i] = 2 * lx[i] + 1;
  clab_webpt(lx, ly, nullptr, 5, nullptr, 1e-15, 100000, o);
  tk_true(T, "web.line", o[0] == 5 && o[1] == 11 && o[5] == 2);
  // a zero weight removes its point
  {
    double zx[4] = {0, 4, 1, 50}, zy[4] = {0, 0, 3, 50}, zw[4] = {1, 1, 1, 0},
           o2[6];
    clab_webpt(zx, zy, zw, 4, nullptr, 1e-14, 100000, o);
    clab_webpt(zx, zy, nullptr, 3, nullptr, 1e-14, 100000, o2);
    tk_eq(T, "web.zerow", o[0], o2[0], 1e-9);
    tk_eq(T, "web.zerowy", o[1], o2[1], 1e-9);
  }
  // a point holding more than half the weight is the minimiser (exactly
  // half would make a whole segment minimal)
  double mw[5] = {1, 1, 1, 1, 5};
  clab_webpt(lx, ly, mw, 5, nullptr, 1e-15, 100000, o);
  tk_true(T, "web.majx", o[0] == 20 && o[1] == 41 && o[5] == 2);
  // coincident points add their weights
  double dx[6] = {0, 0, 0, 3, -2, 1}, dy[6] = {0, 0, 0, 1, 4, -3};
  clab_webpt(dx, dy, nullptr, 6, nullptr, 1e-15, 100000, o);
  tk_true(T, "web.dup", o[0] == 0 && o[1] == 0 && o[5] == 2);
  // parity with pcaPP 2.0.5 l1median_NLM(cbind(rep(x, w), rep(y, w)),
  // maxit = 1000, tol = 1e-14)$par on set A (dev/ref/cpp/cgr.R)
  setA();
  clab_webpt(ax, ay, aw, NA, nullptr, 1e-14, 100000, o);
  tk_eq(T, "web.pcappx", o[0], 501.63658252627886, 1e-8);
  tk_eq(T, "web.pcappy", o[1], 515.7280073006508, 1e-8);
  tk_true(T, "web.pcappobj", o[2] <= 152782.98631206047 * (1 + 1e-14));
  tk_true(T, "web.pcappcode", o[5] >= 1);
  // random sets: converged, stationary, and equal to a compass search
  cgeo::rng r(5);
  for (int rep = 0; rep < 30; ++rep) {
    const int n = 3 + rep * 5;
    double px[200], py[200], pw[200];
    for (int i = 0; i < n; ++i) {
      px[i] = r.unif() * 10 + (i % 4 == 0 ? 20 : 0);
      py[i] = r.unif() * 10;
      pw[i] = 0.5 + r.unif();
    }
    clab_webpt(px, py, pw, n, nullptr, 1e-13, 100000, o);
    tk_true(T, "web.conv", o[5] >= 1);
    double W = 0;
    for (int i = 0; i < n; ++i) W += pw[i];
    tk_true(T, "web.stat", o[3] <= 1e-8 * W);
    if (rep < 10) {
      double qx = o[0] + 1, qy = o[1] - 1;
      compass(px, py, pw, n, &qx, &qy, 1.0);
      tk_eq(T, "web.compass", o[2], wobj(px, py, pw, n, qx, qy), 1e-12);
      tk_eq(T, "web.objcheck", o[2], wobj(px, py, pw, n, o[0], o[1]), 1e-13);
    }
    // monotone: one step at a time from a poor start
    if (rep % 5 == 0) {
      double s[2] = {px[0] + 3, py[0] - 3}, prev = 1e300;
      bool mono = true;
      for (int it = 0; it < 100; ++it) {
        clab_webpt(px, py, pw, n, s, 1e-300, 1, o);
        if (o[2] > prev * (1 + 1e-15)) mono = false;
        prev = o[2];
        s[0] = o[0];
        s[1] = o[1];
      }
      tk_true(T, "web.mono", mono);
    }
    // rotation, scaling, translation, permutation
    if (rep % 3 == 0) {
      double a[6], qx[200], qy[200], qw[200];
      clab_webpt(px, py, pw, n, nullptr, 1e-14, 100000, a);
      const double c = std::cos(1.1), s = std::sin(1.1);
      for (int i = 0; i < n; ++i) {
        qx[n - 1 - i] = 1e3 * (c * px[i] - s * py[i]) + 1e6;
        qy[n - 1 - i] = 1e3 * (s * px[i] + c * py[i]) - 2e6;
        qw[n - 1 - i] = pw[i];
      }
      clab_webpt(qx, qy, qw, n, nullptr, 1e-14, 100000, o);
      tk_eq(T, "web.invx", o[0], 1e3 * (c * a[0] - s * a[1]) + 1e6, 1e-9);
      tk_eq(T, "web.invy", o[1], 1e3 * (s * a[0] + c * a[1]) - 2e6, 1e-9);
    }
  }
  // guards
  double one = 1;
  tk_true(T, "web.n1", clab_webpt(&one, &one, nullptr, 1, nullptr, 1e-9, 10, o) ==
                       CLAB_OK &&
                   o[0] == 1 && o[1] == 1 && o[5] == 2);
  tk_true(T, "web.tol", clab_webpt(sx, sy, nullptr, 4, nullptr, 0, 10, o) ==
                    CLAB_EDOM);
  tk_true(T, "web.maxit", clab_webpt(sx, sy, nullptr, 4, nullptr, 1e-9, -1, o) ==
                      CLAB_EARG);
  double bs[2] = {std::nan(""), 0};
  tk_true(T, "web.start", clab_webpt(sx, sy, nullptr, 4, bs, 1e-9, 10, o) ==
                      CLAB_EDOM);
}

void t_sdell(tk &T, int big) {
  double o[6];
  // points on an ellipse with semi-axes 3, 1 rotated by 0.4, centre (7, -2),
  // 16 equally spaced parameters: the covariance (k = 0) is
  // R diag(a^2 / 2, b^2 / 2) R^T
  double x[16], y[16];
  for (int i = 0; i < 16; ++i) {
    const double t = 2 * PI * i / 16, u = 3 * std::cos(t), v = std::sin(t);
    x[i] = 7 + std::cos(0.4) * u - std::sin(0.4) * v;
    y[i] = -2 + std::sin(0.4) * u + std::cos(0.4) * v;
  }
  tk_true(T, "sde.ok", clab_sdell(x, y, nullptr, 16, 0, 1.0, o) == CLAB_OK);
  tk_eq(T, "sde.mx", o[0], 7, 1e-14);
  tk_eq(T, "sde.a", o[2], 3 / std::sqrt(2.0), 1e-13);
  tk_eq(T, "sde.b", o[3], 1 / std::sqrt(2.0), 1e-13);
  tk_eq(T, "sde.t", o[4], 0.4, 1e-13);
  tk_eq(T, "sde.area", o[5], PI * 1.5, 1e-13);
  // the rotated-axes formulas (Ebdon 1988) with divisor N - 2 against
  // k = 2: the same two standard deviations
  setA();
  {
    double mx = 0, my = 0;
    for (int i = 0; i < NA; ++i) {
      mx += ax[i];
      my += ay[i];
    }
    mx /= NA;
    my /= NA;
    double sxx = 0, syy = 0, sxy = 0;
    for (int i = 0; i < NA; ++i) {
      sxx += (ax[i] - mx) * (ax[i] - mx);
      syy += (ay[i] - my) * (ay[i] - my);
      sxy += (ax[i] - mx) * (ay[i] - my);
    }
    const double th = std::atan(((sxx - syy) + std::sqrt((sxx - syy) * (sxx - syy) +
                                                       4 * sxy * sxy)) /
                                (2 * sxy));
    double a = 0, b = 0;
    for (int i = 0; i < NA; ++i) {
      const double u = ax[i] - mx, v = ay[i] - my;
      a += std::pow(u * std::cos(th) - v * std::sin(th), 2);
      b += std::pow(u * std::sin(th) + v * std::cos(th), 2);
    }
    a = std::sqrt(a / (NA - 2));
    b = std::sqrt(b / (NA - 2));
    clab_sdell(ax, ay, nullptr, NA, 2, 1.0, o);
    tk_eq(T, "sde.ebdonA", o[2], a > b ? a : b, 1e-12);
    tk_eq(T, "sde.ebdonB", o[3], a > b ? b : a, 1e-12);
  }
  // coverage under a bivariate normal: c = sqrt(-2 log(1 - p))
  {
    cgeo::rng r(3);
    const int n = big ? 20000 : 4000;
    static double gx[20000], gy[20000];
    for (int i = 0; i < n; ++i) {
      const double z1 = gauss(r), z2 = gauss(r);
      gx[i] = 2 * z1 + 100;
      gy[i] = 0.75 * z1 + std::sqrt(1 - 0.5625) * z2 - 40;
    }
    const double ps[3] = {0.5, 0.9, 0.39346934028736658};
    for (double p : ps) {
      clab_sdell(gx, gy, nullptr, n, 1, std::sqrt(-2 * std::log(1 - p)), o);
      int in = 0;
      const double c = std::cos(o[4]), s = std::sin(o[4]);
      for (int i = 0; i < n; ++i) {
        const double u = gx[i] - o[0], v = gy[i] - o[1];
        const double a = (c * u + s * v) / o[2], b = (-s * u + c * v) / o[3];
        if (a * a + b * b <= 1) ++in;
      }
      const double f = static_cast<double>(in) / n;
      tk_true(T, "sde.cover", std::fabs(f - p) < 4 * std::sqrt(p * (1 - p) / n));
    }
  }
  // degenerate shapes
  double ix[4] = {2, 2, 2, 2}, iy[4] = {5, 5, 5, 5};
  clab_sdell(ix, iy, nullptr, 4, 1, 1, o);
  tk_true(T, "sde.same", o[2] == 0 && o[3] == 0 && o[4] == 0 && o[5] == 0);
  double lx[5] = {0, 1, 2, 3, 4};
  clab_sdell(lx, lx, nullptr, 5, 1, 1, o);
  tk_eq(T, "sde.lineang", o[4], PI / 4, 1e-15);
  tk_true(T, "sde.lineb", o[3] <= 1e-15 * o[2]);
  double cx[8], cy[8];
  for (int i = 0; i < 8; ++i) {
    cx[i] = std::cos(2 * PI * i / 8);
    cy[i] = std::sin(2 * PI * i / 8);
  }
  clab_sdell(cx, cy, nullptr, 8, 0, 1, o);
  tk_eq(T, "sde.circle", o[2], o[3], 1e-14);
  // rotation equivariance of the angle
  {
    double qx[NA], qy[NA], o2[6];
    clab_sdell(ax, ay, nullptr, NA, 1, 1, o);
    for (int i = 0; i < NA; ++i) {
      qx[i] = std::cos(0.3) * ax[i] - std::sin(0.3) * ay[i];
      qy[i] = std::sin(0.3) * ax[i] + std::cos(0.3) * ay[i];
    }
    clab_sdell(qx, qy, nullptr, NA, 1, 1, o2);
    tk_eq(T, "sde.rotA", o2[2], o[2], 1e-12);
    tk_true(T, "sde.rotT", std::fabs(std::sin(o2[4] - o[4] - 0.3)) < 1e-10);
  }
  // guards
  tk_true(T, "sde.c", clab_sdell(lx, lx, nullptr, 5, 1, 0, o) == CLAB_EDOM);
  tk_true(T, "sde.k", clab_sdell(lx, lx, nullptr, 5, -1, 1, o) == CLAB_EARG);
  tk_true(T, "sde.dof", clab_sdell(lx, lx, nullptr, 2, 2, 1, o) == CLAB_EDOM);
}

void t_cvhull(tk &T, int big) {
  static int idx[2000], hull[4000];
  int nh = 0;
  double o[2];
  // set A against grDevices::chull (R 4.5.0) in dev/ref/cpp/cgr.R: R lists
  // the vertices clockwise. Point 109 is collinear with 187 and 5 in the
  // rationals the set is built from, but not in their doubles: exact
  // rational arithmetic on the doubles (gmp::as.bigq, dev/ref/cpp/cgr.R)
  // puts it 5.4e-12 outside the edge from 187 to 5, so it is a vertex of
  // the hull of the input as given. R's chull drops it; the exact
  // predicate keeps it. Expected: R's list reversed, with 109 inserted.
  setA();
  const int hullA[11] = {99, 60, 8, 0, 145, 5, 109, 187, 192, 132, 33};
  tk_true(T, "hul.ok", clab_cvhull(ax, ay, NA, idx, hull, &nh, o) == CLAB_OK);
  tk_true(T, "hul.nA", nh == 11);
  if (nh == 11) {
    int st = -1;
    for (int k = 0; k < 11; ++k)
      if (hullA[k] == hull[0]) st = k;
    bool same = st >= 0;
    for (int k = 0; k < 11 && same; ++k)
      same = hull[k] == hullA[((st - k) % 11 + 11) % 11];
    tk_true(T, "hul.chull", same);
  }
  // 10 x 10 lattice, every point twice: four corners; R's chull gives the
  // same corners (hullBxy in dev/ref/cpp/cgr.R)
  double bx[200], by[200];
  for (int j = 0; j < 200; ++j) {
    bx[j] = (j % 100) % 10;
    by[j] = (j % 100) / 10;
  }
  clab_cvhull(bx, by, 200, idx, hull, &nh, o);
  tk_true(T, "hul.latn", nh == 4);
  tk_true(T, "hul.latv", hull[0] == 0 && hull[1] == 9 && hull[2] == 99 &&
                     hull[3] == 90);
  tk_eq(T, "hul.latA", o[0], 81, 0);
  tk_eq(T, "hul.latP", o[1], 36, 0);
  // regular heptagon of radius 2
  double px[7], py[7];
  for (int i = 0; i < 7; ++i) {
    px[i] = 2 * std::cos(2 * PI * i / 7) + 1e4;
    py[i] = 2 * std::sin(2 * PI * i / 7) - 5e3;
  }
  clab_cvhull(px, py, 7, idx, hull, &nh, o);
  tk_true(T, "hul.hep", nh == 7);
  tk_eq(T, "hul.hepA", o[0], 3.5 * 4 * std::sin(2 * PI / 7), 1e-11);
  tk_eq(T, "hul.hepP", o[1], 14 * 2 * std::sin(PI / 7), 1e-12);
  // collinear points: two ends, area 0, perimeter twice the length
  double lx[6] = {3, 1, 2, 5, 4, 0}, ly[6];
  for (int i = 0; i < 6; ++i) ly[i] = 0.5 * lx[i];
  clab_cvhull(lx, ly, 6, idx, hull, &nh, o);
  tk_true(T, "hul.line", nh == 2 && hull[0] == 5 && hull[1] == 3 && o[0] == 0);
  tk_eq(T, "hul.lineP", o[1], 2 * std::hypot(5, 2.5), 1e-15);
  double one[2] = {4, 4};
  clab_cvhull(one, one + 1, 1, idx, hull, &nh, o);
  tk_true(T, "hul.one", nh == 1 && o[0] == 0 && o[1] == 0);
  double same[3] = {2, 2, 2};
  clab_cvhull(same, same, 3, idx, hull, &nh, o);
  tk_true(T, "hul.same", nh == 1 && hull[0] == 0);
  // the exact predicate on Kettner et al.'s near-degenerate configuration
  // (doi:10.1016/j.comgeo.2007.06.003): q = (12, 12), r = (24, 24) and the
  // 256 x 256 grid 0.5 + 2^-53 (i, j); integer arithmetic gives the truth
  {
    const double u = 0x1p-53;
    int bad = 0, naive = 0;
    for (int i = 0; i < 256; ++i)
      for (int j = 0; j < 256; ++j) {
        const double qx = 0.5 + i * u, qy = 0.5 + j * u;
        const int t = iorient(qx, qy, 12, 12, 24, 24);
        if (cgeo::orient(qx, qy, 12, 12, 24, 24) != t) ++bad;
        if (norient(qx, qy, 12, 12, 24, 24) != t) ++naive;
      }
    tk_true(T, "orient.exact", bad == 0);
    info(big, "  naive orientation wrong on %d of 65536 Kettner cases\n",
                naive);
    // random triples on the 2^-53 grid near 0.5 with far points
    cgeo::rng r(17);
    bad = 0;
    for (int k = 0; k < (big ? 200000 : 20000); ++k) {
      const double a = 0.5 + r.below(1 << 20) * u, b = 0.5 + r.below(1 << 20) * u;
      const double c = 1 + r.below(31), d = c + r.below(3) * 0.5;
      const double e = 0.5 + r.below(1 << 20) * u, f = 0.5 + r.below(1 << 20) * u;
      if (cgeo::orient(a, b, c, d, e, f) != iorient(a, b, c, d, e, f)) ++bad;
      if (cgeo::orient(c, d, a, b, e, f) != iorient(c, d, a, b, e, f)) ++bad;
    }
    tk_true(T, "orient.random", bad == 0);
    // hull of the Kettner set against gift wrapping with the integer
    // predicate
    static double kx[1026], ky[1026];
    int n = 0;
    kx[n] = 12;
    ky[n++] = 12;
    kx[n] = 24;
    ky[n++] = 24;
    for (int i = 0; i < 32; ++i)
      for (int j = 0; j < 32; ++j) {
        kx[n] = 0.5 + i * u * 7;
        ky[n++] = 0.5 + j * u * 5;
      }
    static int ref[1026];
    const int nr = giftwrap(kx, ky, n, iorient, ref);
    clab_cvhull(kx, ky, n, idx, hull, &nh, o);
    bool eq = nh == nr;
    for (int k = 0; k < nh && eq; ++k) eq = hull[k] == ref[k];
    tk_true(T, "hul.kettner", eq);
  }
  // random points: monotone chain against gift wrapping, and permutation
  {
    cgeo::rng r(23);
    for (int rep = 0; rep < 20; ++rep) {
      const int n = 50 + 40 * rep;
      static double qx[2000], qy[2000], sx_[2000], sy_[2000];
      for (int i = 0; i < n; ++i) {
        // integer coordinates: duplicates and collinear triples occur
        qx[i] = r.below(60);
        qy[i] = r.below(60);
      }
      // gift wrap needs distinct points: deduplicate by brute force
      int m = 0;
      static int keep[2000];
      for (int i = 0; i < n; ++i) {
        bool dup = false;
        for (int j = 0; j < i && !dup; ++j) dup = qx[j] == qx[i] && qy[j] == qy[i];
        if (!dup) keep[m++] = i;
      }
      for (int k = 0; k < m; ++k) {
        sx_[k] = qx[keep[k]];
        sy_[k] = qy[keep[k]];
      }
      static int ref[2000];
      const int nr = giftwrap(sx_, sy_, m, iorient, ref);
      clab_cvhull(qx, qy, n, idx, hull, &nh, o);
      bool eq = nh == nr;
      for (int k = 0; k < nh && eq; ++k)
        eq = qx[hull[k]] == sx_[ref[k]] && qy[hull[k]] == sy_[ref[k]];
      tk_true(T, "hul.brute", eq);
      // reversed input order: the same vertex positions
      static double tx[2000], ty[2000];
      for (int i = 0; i < n; ++i) {
        tx[i] = qx[n - 1 - i];
        ty[i] = qy[n - 1 - i];
      }
      static int h2[4000];
      int nh2 = 0;
      double o2[2];
      clab_cvhull(tx, ty, n, idx, h2, &nh2, o2);
      bool pe = nh2 == nh && o2[0] == o[0];
      for (int k = 0; k < nh && pe; ++k)
        pe = tx[h2[k]] == qx[hull[k]] && ty[h2[k]] == qy[hull[k]];
      tk_true(T, "hul.perm", pe);
      // shoelace against the brute-force vertex list
      double a = 0;
      for (int k = 0; k < nr; ++k) {
        const int p = ref[k], q = ref[(k + 1) % nr];
        a += sx_[p] * sy_[q] - sx_[q] * sy_[p];
      }
      tk_eq(T, "hul.area", o[0], nr >= 3 ? a / 2 : 0, 1e-14);
    }
  }
  // guards
  double huge = 1e200, nanv = std::nan(""), z = 0;
  tk_true(T, "hul.big", clab_cvhull(&huge, &z, 1, idx, hull, &nh, o) == CLAB_EDOM);
  tk_true(T, "hul.nan", clab_cvhull(&nanv, &z, 1, idx, hull, &nh, o) == CLAB_EDOM);
  tk_true(T, "hul.n0", clab_cvhull(&z, &z, 0, idx, hull, &nh, o) == CLAB_EARG);
}

void famtest(tk &T, int big) {
  t_cgstat(T, big);
  t_cgmedn(T, big);
  t_webpt(T, big);
  t_sdell(T, big);
  t_cvhull(T, big);
}

#ifdef CLAB_MAIN
double now() {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

// a heap buffer for the timings (kept out of the program's static size)
template <class T>
struct hb {
  T *p;
  explicit hb(long n) : p(static_cast<T *>(std::malloc(sizeof(T) * n))) {
    if (p == nullptr) std::abort();
  }
  ~hb() { std::free(p); }
  hb(const hb &) = delete;
  hb &operator=(const hb &) = delete;
};

double med5(double *v) { return tk_q7(v, 5, 0.5); }

void fambench() {
  const int ns[3] = {1000, 100000, 1000000};
  hb<double> x_(1000000), y_(1000000);
  double *x = x_.p, *y = y_.p;
  hb<int> idx_(1000000), hull_(2000000);
  int *idx = idx_.p, *hull = hull_.p;
  cgeo::rng r(1);
  for (int i = 0; i < 1000000; ++i) {
    x[i] = r.unif() * 1e4 + 6e5;
    y[i] = r.unif() * 1e4 + 4.8e6;
  }
  for (int n : ns) {
    double t[4][5], o[10];
    int nh;
    for (int k = 0; k < 5; ++k) {
      double a = now();
      clab_cgstat(x, y, nullptr, n, 1, o);
      t[0][k] = now() - a;
      a = now();
      clab_sdell(x, y, nullptr, n, 1, 1, o);
      t[1][k] = now() - a;
      a = now();
      clab_webpt(x, y, nullptr, n, nullptr, 1e-10, 1000, o);
      t[2][k] = now() - a;
      a = now();
      clab_cvhull(x, y, n, idx, hull, &nh, o);
      t[3][k] = now() - a;
    }
    std::printf("bench cgr n = %d: cgstat %.5f sdell %.5f webpt %.5f "
                "cvhull %.5f s (medians of 5)\n",
                n, med5(t[0]), med5(t[1]), med5(t[2]),
                med5(t[3]));
  }
}

#endif

}  // namespace

void t_cgr(tk &T, int big) { famtest(T, big); }
#ifdef CLAB_MAIN
void b_cgr() { fambench(); }
#endif
