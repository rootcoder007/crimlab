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
              double *st) {
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
      // the cluster-robust Wald statistic of the excluded instruments, over q
      std::vector<double> vf(sz(q + px, q + px));
      g = pev::crv(fz.data(), r, w, n, q + px, xi, cl, ncl, cl2, ncl2, c12.data(), n12, 1.0,
                   vf.data(), sc);
      for (int a = 0; a < q; ++a) {
        for (int c = 0; c < q; ++c) xi2[sz(c, q) + static_cast<size_t>(a)] = vf[sz(c, q + px) + static_cast<size_t>(a)];
      }
      if (!pev::chol(xi2, q)) {
        fst[j] = kutil::nan();
      } else {
        std::vector<double> iv(sz(q, q));
        pev::cinv(xi2, q, iv.data());
        double wd = 0.0;
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) wd += b[a] * iv[sz(c, q) + static_cast<size_t>(a)] * b[c];
        }
        fst[j] = wd / q;
      }
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
                int maxit, double eps, double tol, double *beta, double *fst, double *draw,
                double *st) {
  const int p = pe + px + pe;
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
                   nullptr, 0, tol, b1.data(), nullptr, fs, fres.data(), nullptr, s1.data());
    if (rc != CLAB_OK) return rc;
    std::vector<double> xp(sz(m, p));
    for (size_t a = 0; a < sz(m, pe); ++a) xp[a] = e[a];
    for (size_t a = 0; a < sz(m, px); ++a) xp[sz(m, pe) + a] = x2[a];
    for (size_t a = 0; a < sz(m, pe); ++a) xp[sz(m, pe + px) + a] = fres[a];
    std::vector<double> wk(sm * (2 * static_cast<size_t>(p) + 8) + 3 * sz(p, p) +
                           static_cast<size_t>(p) + 4 * (static_cast<size_t>(n1) + static_cast<size_t>(n2)));
    std::vector<int> iw(1);
    return clab_ppmlfe(yy.data(), xp.data(), m, p, off ? oo.data() : nullptr,
                       fw ? ww.data() : nullptr, a1.data(), n1, f2 ? a2.data() : nullptr, n2,
                       nullptr, 0, nullptr, 0, maxit, eps, tol, bt, nullptr, s, wk.data(),
                       iw.data());
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
  std::vector<double> fsb(static_cast<size_t>(pe)), sb(6), bb(static_cast<size_t>(p));
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
                 int base, double *out) {
  if (n < 1 || !(width > 0) || !(width < 360) || !std::isfinite(start) || na < 1) return CLAB_EARG;
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
    int b = static_cast<int>(std::floor(a / width));
    if (b >= nb) b = nb - 1;  // a rounding at the top edge
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

}  // extern "C"
