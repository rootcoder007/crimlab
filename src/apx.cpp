// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, exposure and crime: two-stage least squares with absorbed
// fixed effects (clab_ivfe), the Poisson control function with a cluster
// bootstrap (clab_ppmlcf, on clab_ppmlfe of cnt.cpp), the lookup of the
// same unit on another day (clab_dshift), wind-direction instruments
// (clab_windbin) and the upwind and downwind classification of areas
// around a road or a point source (clab_updwn). The linear algebra is that
// of the policy-evaluation families (pev.h).
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"

namespace {

using std::size_t;

size_t sz(int a, int b) { return static_cast<size_t>(a) * static_cast<size_t>(b); }

constexpr double pi = 3.14159265358979323846;

// The one-effect form of pev::within: its period pass alone.
int demean(double *v, int n, int k, const int *f1, int n1, const int *f2, int n2,
           const double *w, double tol, double *sc) {
  return pev::within(v, n, k, f1, n1, f2 ? f2 : f1, f2 ? n2 : n1, f2 ? 1 : 0, w, tol, 10000, sc);
}

// The scratch of clab_ivfe for n rows, k columns in all, the largest
// regression m columns and effect counts n1, n2.
size_t ivwork(int n, int k, int m, int n1, int n2) {
  const size_t a = 4 * (static_cast<size_t>(n1) + static_cast<size_t>(n2));
  const size_t b = sz(n, m) + 2 * sz(m, m) + static_cast<size_t>(m);
  return sz(n, k) + 3 * static_cast<size_t>(n) + 2 * sz(m, m) + 2 * static_cast<size_t>(m) +
         (a > b ? a : b);
}

// The angle in degrees, 0 to 180, between the directions (ax, ay) and (bx, by).
double ang(double ax, double ay, double bx, double by) {
  const double c = (ax * bx + ay * by) / (std::hypot(ax, ay) * std::hypot(bx, by));
  return std::acos(c > 1 ? 1 : (c < -1 ? -1 : c)) * 180 / pi;
}

}  // namespace

extern "C" {

int clab_ivfe(const double *y, const double *xe, int pe, const double *xx, int px,
              const double *zi, int q, int n, const double *w, const int *f1, int n1,
              const int *f2, int n2, const int *cl, int ncl, const int *cl2, int ncl2,
              double tol, double *beta, double *vc, double *fst, double *fres, double *res,
              double *sco, double *brd, double *st) {
  const int p = pe + px, k = 1 + pe + px + q;
  if (n < 2 || pe < 0 || px < 0 || q < pe || p < 1 || n1 < 1 || (f2 && n2 < 1) ||
      !(tol > 0) || (cl && ncl < 2) || (cl2 && (!cl || ncl2 < 2))) {
    return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if (f1[i] < 0 || f1[i] >= n1 || (f2 && (f2[i] < 0 || f2[i] >= n2))) return CLAB_EARG;
    if ((cl && (cl[i] < 0 || cl[i] >= ncl)) || (cl2 && (cl2[i] < 0 || cl2[i] >= ncl2))) {
      return CLAB_EARG;
    }
    if (w && !(w[i] > 0 && std::isfinite(w[i]))) return CLAB_EDOM;
  }
  const size_t sn = static_cast<size_t>(n);
  const int m = q + px > p ? q + px : p;
  std::vector<double> work(ivwork(n, k, m, n1, n2));
  // v: y, xe, xx, zi, each demeaned
  double *v = work.data(), *r = v + sz(n, k), *u = r + sn;
  double *xi = u + sn, *xi2 = xi + sz(m, m), *b = xi2 + sz(m, m);
  double *sc = b + 2 * static_cast<size_t>(m);
  for (size_t i = 0; i < sn; ++i) v[i] = y[i];
  for (size_t a = 0; a < sz(n, pe); ++a) v[sn + a] = xe[a];
  for (size_t a = 0; a < sz(n, px); ++a) v[sn * (1 + pe) + a] = xx[a];
  for (size_t a = 0; a < sz(n, q); ++a) v[sn * (1 + pe + px) + a] = zi[a];
  const int sw = demean(v, n, k, f1, n1, f2, n2, w, tol, sc);
  if (sw < 0) return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
  double *ye = v, *ee = v + sn, *ex = v + sn * (1 + pe);
  std::vector<int> c12, pw;
  int n12 = 0;
  if (cl && cl2) {
    c12.resize(sn);
    pw.resize(sz(ncl, ncl2));
    n12 = pev::pairs(cl, cl2, n, ncl, ncl2, c12.data(), pw.data());
  }
  int g = 0;
  // H: the second-stage design, the fitted endogenous columns then xx~
  std::vector<double> h(sz(n, p));
  for (size_t a = 0; a < sz(n, px); ++a) h[sz(n, pe) + a] = ex[a];
  // the first stages: each endogenous column on zi~ then xx~
  std::vector<double> fz(sz(n, q + px));
  for (size_t a = 0; a < sz(n, q); ++a) fz[a] = v[sn * (1 + pe + px) + a];
  for (size_t a = 0; a < sz(n, px); ++a) fz[sz(n, q) + a] = ex[a];
  for (int j = 0; j < pe; ++j) {
    const double *ej = ee + sz(n, j);
    if (pev::ols(fz.data(), ej, n, q + px, w, b, xi, r, sc) != CLAB_OK) return CLAB_EDOM;
    for (size_t i = 0; i < sn; ++i) h[sz(n, j) + i] = ej[i] - r[i];
    if (fres) {
      for (size_t i = 0; i < sn; ++i) fres[sz(n, j) + i] = r[i];
    }
    if (cl) {
      // the cluster-robust Wald statistic of the excluded instruments over
      // q; the Kleibergen-Paap rk Wald F, which for one endogenous column
      // is that statistic; and the effective F of Montiel Olea and
      // Pflueger (2013), bz' Q bz / tr(Vzz Q), Q the instruments' cross
      // products after the controls and the effects
      std::vector<double> vf(sz(q + px, q + px)), bz(b, b + q);
      g = pev::crv(fz.data(), r, w, n, q + px, xi, cl, ncl, cl2, ncl2, c12.data(), n12, 1.0,
                   vf.data(), sc);
      std::vector<double> vz(sz(q, q));
      for (int a = 0; a < q; ++a) {
        for (int c = 0; c < q; ++c) {
          vz[sz(c, q) + static_cast<size_t>(a)] = vf[sz(c, q + px) + static_cast<size_t>(a)];
          xi2[sz(c, q) + static_cast<size_t>(a)] = vz[sz(c, q) + static_cast<size_t>(a)];
        }
      }
      double wald = kutil::nan();
      if (pev::chol(xi2, q)) {
        std::vector<double> iv(sz(q, q));
        pev::cinv(xi2, q, iv.data());
        double wd = 0.0;
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) wd += bz[a] * iv[sz(c, q) + static_cast<size_t>(a)] * bz[c];
        }
        wald = wd / q;
      }
      fst[j] = wald;
      fst[pe + j] = pe == 1 ? wald : kutil::nan();
      double eff = kutil::nan();
      if (pe == 1) {
        // the instruments less their fit on the controls (all demeaned)
        std::vector<double> rz(sz(n, q));
        for (int a = 0; a < q; ++a) {
          const double *za = fz.data() + sz(n, a);
          if (px > 0) {
            if (pev::ols(ex, za, n, px, w, b, xi, rz.data() + sz(n, a), sc) != CLAB_OK) return CLAB_EDOM;
          } else {
            for (size_t i = 0; i < sn; ++i) rz[sz(n, a) + i] = za[i];
          }
        }
        std::vector<double> qq(sz(q, q));
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) {
            kutil::ksum t;
            for (size_t i = 0; i < sn; ++i) t.add((w ? w[i] : 1.0) * rz[sz(n, a) + i] * rz[sz(n, c) + i]);
            qq[sz(c, q) + static_cast<size_t>(a)] = t.val();
          }
        }
        double num = 0.0, den = 0.0;
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) {
            num += bz[a] * qq[sz(c, q) + static_cast<size_t>(a)] * bz[c];
            den += vz[sz(c, q) + static_cast<size_t>(a)] * qq[sz(a, q) + static_cast<size_t>(c)];
          }
        }
        if (den > 0) eff = num / den;
      }
      fst[2 * pe + j] = eff;
    }
  }
  // the second stage on H, the residuals on the actual columns
  if (pev::ols(h.data(), ye, n, p, w, beta, xi, r, sc) != CLAB_OK) return CLAB_EDOM;
  for (size_t i = 0; i < sn; ++i) {
    double f = 0.0;
    for (int j = 0; j < pe; ++j) f += ee[sz(n, j) + i] * beta[j];
    for (int j = 0; j < px; ++j) f += ex[sz(n, j) + i] * beta[pe + j];
    u[i] = ye[i] - f;
    if (res) res[i] = u[i];
  }
  if (brd) {
    for (size_t a = 0; a < sz(p, p); ++a) brd[a] = xi[a];
  }
  if (sco) {
    for (int j = 0; j < p; ++j) {
      for (size_t i = 0; i < sn; ++i) sco[sz(n, j) + i] = (w ? w[i] : 1.0) * u[i] * h[sz(n, j) + i];
    }
  }
  if (cl && vc) {
    g = pev::crv(h.data(), u, w, n, p, xi, cl, ncl, cl2, ncl2, c12.data(), n12, 1.0, vc, sc);
  }
  kutil::ksum rss;
  for (size_t i = 0; i < sn; ++i) rss.add((w ? w[i] : 1.0) * u[i] * u[i]);
  st[0] = g;
  st[1] = sw;
  st[2] = rss.val();
  return CLAB_OK;
}

int clab_ppmlcf(const double *y, const double *xe, int pe, const double *xx, int px,
                const double *zi, int q, int n, const double *off, const int *f1, int n1,
                const int *f2, int n2, const int *cl, int ncl, const int *bm, int nb,
                int how, int maxit, double eps, double tol, double *beta, double *fst,
                double *draw, double *st) {
  // how 0: the control function (the residuals join the regressors); 1:
  // the fitted exposure in place of the exposure
  const int p = how == 1 ? pe + px : pe + px + pe;
  if (how < 0 || how > 1) return CLAB_EARG;
  if (n < 2 || pe < 1 || px < 0 || q < pe || ncl < 2 || nb < 0 || !cl) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (cl[i] < 0 || cl[i] >= ncl) return CLAB_EARG;
  }
  // one fit on the rows with weights fw (NULL: all rows, weight 1); the
  // first stage's residuals join the regressors of the Poisson fit
  auto one = [&](const std::vector<int> &rows, const double *fw, const int *c, int nc, double *bt,
                 double *fs, double *s) -> int {
    const int m = static_cast<int>(rows.size());
    const size_t sm = static_cast<size_t>(m);
    std::vector<double> yy(sm), ww, oo, e(sz(m, pe)), x2(sz(m, px)), z2(sz(m, q));
    std::vector<int> a1(sm), a2(f2 ? sm : 0), cc(sm);
    for (size_t i = 0; i < sm; ++i) {
      const int r = rows[i];
      yy[i] = y[r];
      a1[i] = f1[r];
      if (f2) a2[i] = f2[r];
      cc[i] = c ? c[r] : 0;
      for (int j = 0; j < pe; ++j) e[sz(m, j) + i] = xe[sz(n, j) + static_cast<size_t>(r)];
      for (int j = 0; j < px; ++j) x2[sz(m, j) + i] = xx[sz(n, j) + static_cast<size_t>(r)];
      for (int j = 0; j < q; ++j) z2[sz(m, j) + i] = zi[sz(n, j) + static_cast<size_t>(r)];
    }
    if (fw) {
      ww.resize(sm);
      for (size_t i = 0; i < sm; ++i) ww[i] = fw[rows[i]];
    }
    if (off) {
      oo.resize(sm);
      for (size_t i = 0; i < sm; ++i) oo[i] = off[rows[i]];
    }
    std::vector<double> fres(sz(m, pe)), b1(static_cast<size_t>(pe + px)), s1(3);
    const int rc = clab_ivfe(yy.data(), e.data(), pe, x2.data(), px, z2.data(), q, m, fw ? ww.data() : nullptr,
                   a1.data(), n1, f2 ? a2.data() : nullptr, n2, c ? cc.data() : nullptr, nc,
                   nullptr, 0, tol, b1.data(), nullptr, fs, fres.data(), nullptr, nullptr, nullptr,
                   s1.data());
    if (rc != CLAB_OK) return rc;
    std::vector<double> xp(sz(m, p));
    for (size_t a = 0; a < sz(m, pe); ++a) xp[a] = how == 1 ? e[a] - fres[a] : e[a];
    for (size_t a = 0; a < sz(m, px); ++a) xp[sz(m, pe) + a] = x2[a];
    if (how == 0) {
      for (size_t a = 0; a < sz(m, pe); ++a) xp[sz(m, pe + px) + a] = fres[a];
    }
    std::vector<double> wk(sm * (2 * static_cast<size_t>(p) + 8) + 3 * sz(p, p) +
                           static_cast<size_t>(p) + 4 * (static_cast<size_t>(n1) + static_cast<size_t>(n2)));
    std::vector<int> iw(1);
    return clab_ppmlfe(yy.data(), xp.data(), m, p, off ? oo.data() : nullptr,
                       fw ? ww.data() : nullptr, a1.data(), n1, f2 ? a2.data() : nullptr, n2,
                       nullptr, 0, nullptr, 0, maxit, eps, tol, bt, nullptr, nullptr, nullptr, s,
                       wk.data(), iw.data());
  };
  std::vector<int> all(static_cast<size_t>(n));
  for (int i = 0; i < n; ++i) all[static_cast<size_t>(i)] = i;
  double s0[6];
  int rc = one(all, nullptr, cl, ncl, beta, fst, s0);
  if (rc != CLAB_OK) return rc;
  st[0] = s0[0];
  st[1] = s0[1];
  int bad = 0;
  std::vector<double> fw(static_cast<size_t>(n)), sep(static_cast<size_t>(n1 + n2));
  std::vector<double> fsb(3 * static_cast<size_t>(pe)), sb(6), bb(static_cast<size_t>(p));
  for (int d = 0; d < nb; ++d) {
    if (kutil::stop(static_cast<double>(n) * p * 20)) return CLAB_EINTR;
    // the rows of the clusters drawn, each weighted by how often it was drawn
    std::vector<int> rows;
    for (int i = 0; i < n; ++i) {
      fw[static_cast<size_t>(i)] = bm[sz(d, ncl) + static_cast<size_t>(cl[i])];
      if (fw[static_cast<size_t>(i)] > 0) rows.push_back(i);
    }
    // a group of zeros in this draw is dropped, as in the full sample
    std::vector<double> yr(rows.size());
    std::vector<int> r1(rows.size()), r2(f2 ? rows.size() : 0), kp(rows.size());
    for (size_t a = 0; a < rows.size(); ++a) {
      yr[a] = y[rows[a]];
      r1[a] = f1[rows[a]];
      if (f2) r2[a] = f2[rows[a]];
    }
    std::vector<int> kept;
    if (!rows.empty() &&
        clab_fesep(yr.data(), static_cast<int>(rows.size()), r1.data(), n1, f2 ? r2.data() : nullptr,
                   n2, kp.data(), sep.data()) >= 0) {
      for (size_t a = 0; a < rows.size(); ++a) {
        if (kp[a]) kept.push_back(rows[a]);
      }
    }
    rc = kept.size() < 2 ? CLAB_EDOM
                         : one(kept, fw.data(), nullptr, 0, bb.data(), fsb.data(), sb.data());
    if (rc == CLAB_EINTR) return rc;
    for (int j = 0; j < p; ++j) {
      draw[sz(j, nb) + static_cast<size_t>(d)] = rc == CLAB_OK ? bb[static_cast<size_t>(j)] : kutil::nan();
    }
    bad += rc != CLAB_OK;
  }
  st[2] = bad;
  return CLAB_OK;
}

int clab_windbin(const double *dir, int n, double width, double start, const int *area, int na,
                 int base, int edge, double *out) {
  if (n < 1 || !(width > 0) || !(width < 360) || !std::isfinite(start) || na < 1 || edge < 0 ||
      edge > 1) {
    return CLAB_EARG;
  }
  const double nbd = 360 / width;
  const int nb = static_cast<int>(std::nearbyint(nbd));
  if (std::fabs(nbd - nb) > 1e-9 || nb < 2 || base < 0 || base >= nb) return CLAB_EARG;
  const int k = nb - 1;
  for (size_t a = 0; a < sz(n, na * k); ++a) out[a] = 0.0;
  for (int i = 0; i < n; ++i) {
    const int ar = area ? area[i] : 0;
    if (ar < 0 || ar >= na) return CLAB_EARG;
    if (!std::isfinite(dir[i])) {
      for (int c = 0; c < na * k; ++c) out[sz(c, n) + static_cast<size_t>(i)] = kutil::nan();
      continue;
    }
    double a = std::fmod(dir[i] - start, 360.0);
    if (a < 0) a += 360.0;
    // edge 0: [lo, hi); edge 1: (lo, hi], the start itself in the first
    // sector (whole degrees 0-90, 91-180, ... as Bondy, Roth and Sager)
    int b = edge ? (a == 0 ? 0 : static_cast<int>(std::ceil(a / width)) - 1)
                 : static_cast<int>(std::floor(a / width));
    if (b >= nb) b = nb - 1;  // a rounding at the top edge
    if (b < 0) b = 0;
    if (b == base) continue;
    const int c = ar * k + (b < base ? b : b - 1);
    out[sz(c, n) + static_cast<size_t>(i)] = 1.0;
  }
  return CLAB_OK;
}

int clab_dshift(const int *u, const double *t, int n, double k, int *idx, int *ord) {
  if (n < 1 || !std::isfinite(k)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(t[i])) return CLAB_EDOM;
    ord[i] = i;
  }
  auto less = [&](int a, int b) { return u[a] < u[b] || (u[a] == u[b] && t[a] < t[b]); };
  kutil::hsort(ord, n, less);
  for (int q = 1; q < n; ++q) {
    if (u[ord[q]] == u[ord[q - 1]] && t[ord[q]] == t[ord[q - 1]]) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    const double want = t[i] + k;
    int lo = 0, hi = n;  // the first position not below (u[i], want)
    while (lo < hi) {
      const int mid = lo + (hi - lo) / 2;
      const int o = ord[mid];
      if (u[o] < u[i] || (u[o] == u[i] && t[o] < want)) {
        lo = mid + 1;
      } else {
        hi = mid;
      }
    }
    idx[i] = lo < n && u[ord[lo]] == u[i] && t[ord[lo]] == want ? ord[lo] : -1;
  }
  return CLAB_OK;
}

int clab_updwn(const double *ax, const double *ay, const double *wd, int n, const double *sx0,
               const double *sy0, const double *sx1, const double *sy1, int ns, double band,
               double tol, int twin, int *src, int *side, int *dn, double *dst) {
  if (n < 1 || ns < 1 || !(band > 0) || !(tol > 0) || tol > 90 || twin < 0 || twin > 1) {
    return CLAB_EARG;
  }
  for (int s = 0; s < ns; ++s) {
    if (!std::isfinite(sx0[s]) || !std::isfinite(sy0[s])) return CLAB_EDOM;
    const bool line = std::isfinite(sx1[s]) && std::isfinite(sy1[s]);
    if (line && sx1[s] == sx0[s] && sy1[s] == sy0[s]) return CLAB_EDOM;
  }
  for (int i = 0; i < n; ++i) {
    src[i] = -1;
    side[i] = -1;
    dn[i] = -2;
    if (dst) dst[i] = kutil::nan();
    if (!std::isfinite(ax[i]) || !std::isfinite(ay[i]) || !std::isfinite(wd[i])) {
      dn[i] = -5;
      continue;
    }
    // the sources within band: exactly one is needed
    int hit = -1, nh = 0;
    double best = 0.0;
    for (int s = 0; s < ns; ++s) {
      double d;
      if (std::isfinite(sx1[s]) && std::isfinite(sy1[s])) {
        const double dx = sx1[s] - sx0[s], dy = sy1[s] - sy0[s];
        double t = ((ax[i] - sx0[s]) * dx + (ay[i] - sy0[s]) * dy) / (dx * dx + dy * dy);
        t = t < 0 ? 0 : (t > 1 ? 1 : t);
        d = std::hypot(ax[i] - sx0[s] - t * dx, ay[i] - sy0[s] - t * dy);
      } else {
        d = std::hypot(ax[i] - sx0[s], ay[i] - sy0[s]);
      }
      if (d <= band) {
        // twin 1 keeps the nearest (the first of equal distances)
        if (nh == 0 || d < best) {
          hit = s;
          best = d;
        }
        ++nh;
      }
    }
    if (nh == 0) continue;
    if (nh > 1 && twin == 0) {
      dn[i] = -3;
      continue;
    }
    src[i] = hit;
    if (dst) dst[i] = best;
    // where the wind blows to: from direction theta, clockwise from north,
    // x east and y north
    const double th = wd[i] * pi / 180;
    const double tx = -std::sin(th), ty = -std::cos(th);
    double nx, ny;
    if (std::isfinite(sx1[hit]) && std::isfinite(sy1[hit])) {
      const double dx = sx1[hit] - sx0[hit], dy = sy1[hit] - sy0[hit];
      const double cr = dx * (ay[i] - sy0[hit]) - dy * (ax[i] - sx0[hit]);
      if (cr == 0) {
        dn[i] = -4;
        continue;
      }
      // the normal of the line towards the area's side: left is (-dy, dx)
      side[i] = cr > 0 ? 1 : 0;
      nx = cr > 0 ? -dy : dy;
      ny = cr > 0 ? dx : -dx;
    } else {
      nx = ax[i] - sx0[hit];
      ny = ay[i] - sy0[hit];
      if (nx == 0 && ny == 0) {
        dn[i] = -4;
        continue;
      }
    }
    const double a = ang(tx, ty, nx, ny);
    dn[i] = a <= tol ? 1 : (a >= 180 - tol ? 0 : -1);
  }
  return CLAB_OK;
}

int clab_vcx(const double *sco, int n, int p, const double *brd, int kind, const int *cl, int ncl,
             const int *cl2, int ncl2, const double *t, int lag, double fac, double *vc) {
  if (n < 2 || p < 1 || kind < 0 || kind > 3 || !(fac >= 0)) return CLAB_EARG;
  if ((kind == 0 || kind == 1 || kind == 3) && (!cl || ncl < 1)) return CLAB_EARG;
  if (kind == 1 && (!cl2 || ncl2 < 1)) return CLAB_EARG;
  if (kind == 3 && (!t || lag < 0)) return CLAB_EARG;
  const size_t sn = static_cast<size_t>(n), sp = static_cast<size_t>(p);
  std::vector<double> m(sp * sp, 0.0);
  // m += sg * sum over the groups g of s_g s_g'
  auto grp = [&](const int *c, int nc, double sg) {
    std::vector<double> s(static_cast<size_t>(nc) * sp, 0.0);
    for (size_t j = 0; j < sp; ++j) {
      for (size_t i = 0; i < sn; ++i) s[static_cast<size_t>(c[i]) * sp + j] += sco[j * sn + i];
    }
    for (int g = 0; g < nc; ++g) {
      const double *v = s.data() + static_cast<size_t>(g) * sp;
      for (size_t j = 0; j < sp; ++j) {
        for (size_t k = 0; k < sp; ++k) m[k * sp + j] += sg * v[j] * v[k];
      }
    }
  };
  for (int i = 0; i < n; ++i) {
    if ((cl && (cl[i] < 0 || cl[i] >= ncl)) || (kind == 1 && (cl2[i] < 0 || cl2[i] >= ncl2))) {
      return CLAB_EARG;
    }
  }
  if (kind == 0) {
    grp(cl, ncl, 1.0);
  } else if (kind == 1) {
    std::vector<int> c12(sn), pw(sz(ncl, ncl2));
    const int n12 = pev::pairs(cl, cl2, n, ncl, ncl2, c12.data(), pw.data());
    grp(cl, ncl, 1.0);
    grp(cl2, ncl2, 1.0);
    grp(c12.data(), n12, -1.0);
  } else {
    // each row alone (kind 2), then for Newey-West (kind 3) the pairs of
    // rows of one unit whose times differ by 1 to lag, with Bartlett
    // weights 1 - d / (lag + 1)
    for (size_t i = 0; i < sn; ++i) {
      for (size_t j = 0; j < sp; ++j) {
        for (size_t k = 0; k < sp; ++k) m[k * sp + j] += sco[j * sn + i] * sco[k * sn + i];
      }
    }
    if (kind == 3 && lag > 0) {
      std::vector<int> ord(sn);
      for (int i = 0; i < n; ++i) {
        if (!std::isfinite(t[i])) return CLAB_EDOM;
        ord[static_cast<size_t>(i)] = i;
      }
      auto less = [&](int a, int b) { return cl[a] < cl[b] || (cl[a] == cl[b] && t[a] < t[b]); };
      kutil::hsort(ord.data(), n, less);
      for (size_t a = 0; a < sn; ++a) {
        const int i = ord[a];
        for (size_t b = a + 1; b < sn; ++b) {
          const int r = ord[b];
          if (cl[r] != cl[i]) break;
          const double d = t[r] - t[i];
          if (d > lag) break;
          if (!(d > 0)) continue;
          const double wl = 1.0 - d / (lag + 1.0);
          for (size_t j = 0; j < sp; ++j) {
            for (size_t k = 0; k < sp; ++k) {
              m[k * sp + j] += wl * (sco[j * sn + static_cast<size_t>(i)] * sco[k * sn + static_cast<size_t>(r)] +
                                     sco[j * sn + static_cast<size_t>(r)] * sco[k * sn + static_cast<size_t>(i)]);
            }
          }
        }
      }
    }
  }
  // fac brd m brd
  std::vector<double> t1(sp * sp, 0.0);
  for (size_t j = 0; j < sp; ++j) {
    for (size_t k = 0; k < sp; ++k) {
      double a = 0.0;
      for (size_t l = 0; l < sp; ++l) a += brd[l * sp + j] * m[k * sp + l];
      t1[k * sp + j] = a;
    }
  }
  for (size_t j = 0; j < sp; ++j) {
    for (size_t k = 0; k < sp; ++k) {
      double a = 0.0;
      for (size_t l = 0; l < sp; ++l) a += t1[l * sp + j] * brd[k * sp + l];
      vc[k * sp + j] = fac * a;
    }
  }
  return CLAB_OK;
}

}  // extern "C"
