// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The planar geometry declared in cgeo.h, compiled once.
#include "cgeo.h"

#include <cmath>

#include "kutil.h"

namespace cgeo {

namespace {

// Error-free transformations (Dekker 1971, doi:10.1007/BF01397083; Knuth).
// Every product inside two_prod is exact, so a compiler that fuses a
// multiply and an add there cannot change the result.
void two_sum(double a, double b, double *s, double *e) {
  const double x = a + b;
  const double bv = x - a;
  const double av = x - bv;
  *s = x;
  *e = (a - av) + (b - bv);
}
void split(double a, double *hi, double *lo) {
  const double c = 134217729.0 * a;  // 2^27 + 1
  const double ab = c - a;
  *hi = c - ab;
  *lo = a - *hi;
}
void two_prod(double a, double b, double *p, double *e) {
  const double x = a * b;
  double ah, al, bh, bl;
  split(a, &ah, &al);
  split(b, &bh, &bl);
  *p = x;
  *e = ((ah * bh - x) + ah * bl + al * bh) + al * bl;
}
}  // namespace

int orient(double ax, double ay, double bx, double by, double cx,
                  double cy) {
  const double l = (ax - cx) * (by - cy);
  const double r = (ay - cy) * (bx - cx);
  const double det = l - r;
  const double sum = std::fabs(l) + std::fabs(r);
  const double eps = 0x1.0p-53;
  if (sum > 1e-280 && std::fabs(det) > (3.0 + 16.0 * eps) * eps * sum)
    return det > 0.0 ? 1 : -1;
  const double fa[6] = {ax, -ax, -cx, -ay, ay, cy};
  const double fb[6] = {by, cy, by, bx, cx, bx};
  double h[12];
  int m = 0;
  for (int k = 0; k < 6; ++k) {
    double p, e;
    two_prod(fa[k], fb[k], &p, &e);
    const double add[2] = {e, p};
    for (int t = 0; t < 2; ++t) {
      double q = add[t];
      int w = 0;
      for (int j = 0; j < m; ++j) {
        double s, err;
        two_sum(q, h[j], &s, &err);
        q = s;
        if (err != 0.0) h[w++] = err;
      }
      if (q != 0.0) h[w++] = q;
      m = w;
    }
  }
  for (int j = m - 1; j >= 0; --j) {
    if (h[j] > 0.0) return 1;
    if (h[j] < 0.0) return -1;
  }
  return 0;
}

bool moments(const double *x, const double *y, const double *w, int n,
                    double *m) {
  kutil::ksum sw, sx, sy;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (!kutil::nonneg(wi) || !std::isfinite(x[i]) || !std::isfinite(y[i]))
      return false;
    sw.add(wi);
    sx.add(wi * x[i]);
    sy.add(wi * y[i]);
  }
  const double W = sw.val();
  if (!(W > 0.0) || !std::isfinite(W)) return false;
  double mx = sx.val() / W, my = sy.val() / W;
  // one refinement step: the mean of the residuals
  kutil::ksum rx, ry;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    rx.add(wi * (x[i] - mx));
    ry.add(wi * (y[i] - my));
  }
  mx += rx.val() / W;
  my += ry.val() / W;
  kutil::ksum a, b, c;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    const double dx = x[i] - mx, dy = y[i] - my;
    a.add(wi * dx * dx);
    b.add(wi * dy * dy);
    c.add(wi * dx * dy);
  }
  m[0] = W;
  m[1] = mx;
  m[2] = my;
  m[3] = a.val();
  m[4] = b.val();
  m[5] = c.val();
  return true;
}

void eig2(double a, double b, double c, double *l1, double *l2,
                 double *t) {
  const double h = 0.5 * (a - b);
  const double r = std::hypot(h, c);
  const double mid = 0.5 * (a + b);
  *l1 = mid + r;
  const double det = a * b - c * c;
  *l2 = *l1 > 0.0 ? det / *l1 : mid - r;
  // rounding can leave the smaller root of a PSD matrix just below 0
  *l2 = std::fmin(std::fmax(*l2, 0.0), *l1);
  // 0.5 atan2 lies in [-pi / 2, pi / 2]
  double th = (r == 0.0) ? 0.0 : 0.5 * std::atan2(c, h);
  if (th < 0.0) th += pi;
  *t = th;
}

void weber(const double *x, const double *y, const double *w, int n,
                  double ox, double oy, double scale, double tol, int maxit,
                  double *cx, double *cy, double *out) {
  // the objective, the weight sitting at (qx, qy) and the length of the sum
  // of unit vectors from it to the other points
  auto eval = [&](double qx, double qy, double *f, double *eta, double *rx,
                  double *ry, double *s) {
    kutil::ksum S, Rx, Ry, F;
    double e = 0.0;
    for (int i = 0; i < n; ++i) {
      const double wi = w ? w[i] : 1.0;
      if (wi == 0.0) continue;
      const double dx = (x[i] - ox) - qx, dy = (y[i] - oy) - qy;
      const double d = std::hypot(dx, dy);
      if (d == 0.0) {
        e += wi;
        continue;
      }
      S.add(wi / d);
      Rx.add(wi * dx / d);
      Ry.add(wi * dy / d);
      F.add(wi * d);
    }
    *f = F.val();
    *eta = e;
    *rx = Rx.val();
    *ry = Ry.val();
    *s = S.val();
  };
  double px = *cx, py = *cy;
  int it = 0, code = 0;
  double fobj, eta, rx, ry, s;
  for (;;) {
    eval(px, py, &fobj, &eta, &rx, &ry, &s);
    const double r = std::hypot(rx, ry);
    if (s == 0.0 || r <= eta) {
      code = 2;
      break;
    }
    if (it >= maxit) break;
    const double g = eta > 0.0 ? eta / r : 0.0;
    const double stx = (1.0 - g) * rx / s, sty = (1.0 - g) * ry / s;
    px += stx;
    py += sty;
    ++it;
    if (std::hypot(stx, sty) <= tol * scale) {
      code = 1;
      eval(px, py, &fobj, &eta, &rx, &ry, &s);
      break;
    }
  }
  double excess = std::fmax(std::hypot(rx, ry) - eta, 0.0);
  // Iterates that approach a data point which is itself the minimiser
  // converge slowly and never land on it. When the iteration stops short
  // of an exact certificate, the data point nearest the result is tested:
  // if |R| <= eta there, it is a minimiser and is returned (code 2).
  if (code != 2) {
    int q = -1;
    double best = 0.0;
    for (int i = 0; i < n; ++i) {
      const double d = std::hypot((x[i] - ox) - px, (y[i] - oy) - py);
      if ((w ? w[i] : 1.0) > 0.0 && (q < 0 || d < best)) {
        q = i;
        best = d;
      }
    }
    double fq, eq, qx, qy, sq;
    eval(x[q] - ox, y[q] - oy, &fq, &eq, &qx, &qy, &sq);
    if (std::hypot(qx, qy) <= eq) {
      px = x[q] - ox;
      py = y[q] - oy;
      fobj = fq;
      excess = 0.0;
      code = 2;
    }
  }
  *cx = px;
  *cy = py;
  out[0] = fobj;
  out[1] = excess;
  out[2] = it;
  out[3] = code;
}

void knn(const grid &g, const double *x, const double *y, int i, int k,
                double *best) {
  int got = 0;
  const double px = x[i], py = y[i];
  const int ci = g.cx(px), cj = g.cy(py);
  for (int r = 0;; ++r) {
    const int a0 = ci - r, a1 = ci + r, b0 = cj - r, b1 = cj + r;
    for (int b = b0; b <= b1; ++b) {
      if (b < 0 || b >= g.gy) continue;
      const bool edge = (b == b0 || b == b1);
      for (int a = a0; a <= a1; a += (edge ? 1 : a1 - a0 > 0 ? a1 - a0 : 1)) {
        if (a >= 0 && a < g.gx) {
          const int c = b * g.gx + a;
          for (int t = g.start[c]; t < g.start[c + 1]; ++t) {
            const int j = g.ord[t];
            if (j == i) continue;
            const double dx = x[j] - px, dy = y[j] - py;
            const double d = dx * dx + dy * dy;
            if (got < k) {
              int p = got++;
              while (p > 0 && best[p - 1] > d) {
                best[p] = best[p - 1];
                --p;
              }
              best[p] = d;
            } else if (d < best[k - 1]) {
              int p = k - 1;
              while (p > 0 && best[p - 1] > d) {
                best[p] = best[p - 1];
                --p;
              }
              best[p] = d;
            }
          }
        }
        if (a1 == a0) break;
      }
    }
    // every cell of the grid visited
    if (a0 <= 0 && b0 <= 0 && a1 >= g.gx - 1 && b1 >= g.gy - 1) return;
    if (got == k) {
      // the nearest a point outside the block of rings 0..r can be
      double m = px - (g.x0 + a0 * g.wx);
      const double e1 = (g.x0 + (a1 + 1) * g.wx) - px;
      const double e2 = py - (g.y0 + b0 * g.wy);
      const double e3 = (g.y0 + (b1 + 1) * g.wy) - py;
      if (a0 <= 0) m = 1e308;
      double mm = m;
      if (a1 < g.gx - 1 && e1 < mm) mm = e1;
      if (b0 > 0 && e2 < mm) mm = e2;
      if (b1 < g.gy - 1 && e3 < mm) mm = e3;
      if (mm > 0.0 && mm * mm >= best[k - 1]) return;
    }
  }
}

void grid::build(const double *x, const double *y, int n, double ax0, double ax1,
             double ay0, double ay1, double cell, int *iwork) {
    x0 = ax0;
    y0 = ay0;
    const double w = ax1 - ax0, h = ay1 - ay0;
    double c = cell > 0.0 ? cell : 1.0;
    for (;;) {
      const double fx = w > 0.0 ? std::ceil(w / c) : 1.0;
      const double fy = h > 0.0 ? std::ceil(h / c) : 1.0;
      if (fx * fy <= n + 1.0) {
        gx = static_cast<int>(fx < 1.0 ? 1.0 : fx);
        gy = static_cast<int>(fy < 1.0 ? 1.0 : fy);
        break;
      }
      c *= 1.25;
    }
    wx = w > 0.0 ? w / gx : 1.0;
    wy = h > 0.0 ? h / gy : 1.0;
    const int nc = gx * gy;
    start = iwork;
    ord = iwork + nc + 1;
    for (int i = 0; i <= nc; ++i) start[i] = 0;
    for (int i = 0; i < n; ++i) ++start[cy(y[i]) * gx + cx(x[i]) + 1];
    for (int i = 0; i < nc; ++i) start[i + 1] += start[i];
    // place points in index order within each cell; start[c] is used as a
    // cursor and restored afterwards
    for (int i = 0; i < n; ++i) ord[start[cy(y[i]) * gx + cx(x[i])]++] = i;
    for (int i = nc; i > 0; --i) start[i] = start[i - 1];
    start[0] = 0;
  }

}  // namespace cgeo
