// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, exposure and crime: two-stage least squares with absorbed
// fixed effects (clab_ivfe) and its over-identification tests
// (clab_ivoid), the Poisson control function with a cluster bootstrap
// (clab_ppmlcf, on clab_ppmlfe of cnt.cpp), the exponential-mean
// generalised method of moments (clab_expgmm), the lookup of the
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
           const double *w, double tol, double *sc, const int *fx = nullptr,
           const int *nx = nullptr, int mx = 0) {
  return pev::withink(v, n, k, f1, n1, f2, n2, fx, nx, mx, w, tol, 10000, sc);
}

// The effects' codes are in range: f1, f2 (or NULL) and the mx sets fx.
bool fecodes(int n, const int *f1, int n1, const int *f2, int n2, const int *fx, const int *nx,
             int mx) {
  if (n1 < 1 || (f2 && n2 < 1) || mx < 0 || (mx > 0 && (!f2 || !fx || !nx))) return false;
  for (int i = 0; i < n; ++i) {
    if (f1[i] < 0 || f1[i] >= n1 || (f2 && (f2[i] < 0 || f2[i] >= n2))) return false;
    for (int j = 0; j < mx; ++j) {
      const int c = fx[sz(n, j) + static_cast<size_t>(i)];
      if (nx[j] < 1 || c < 0 || c >= nx[j]) return false;
    }
  }
  return true;
}

// n1 + n2 + the sum of nx
int ngrp(int n1, int n2, const int *nx, int mx) {
  int g = n1 + n2;
  for (int j = 0; j < mx; ++j) g += nx[j];
  return g;
}

// The scratch of clab_ivfe for n rows, k columns in all, the largest
// regression m columns and effect counts n1, n2.
size_t ivwork(int n, int k, int m, int n1, int n2) {
  const size_t a = 4 * (static_cast<size_t>(n1) + static_cast<size_t>(n2));
  const size_t b = sz(n, m) + 2 * sz(m, m) + static_cast<size_t>(m);
  return sz(n, k) + 3 * static_cast<size_t>(n) + 2 * sz(m, m) + 2 * static_cast<size_t>(m) +
         (a > b ? a : b);
}

// The inverse of the symmetric positive definite p x p matrix a into ai;
// false when it is singular to working precision. sc: p p doubles.
bool spdinv(const double *a, int p, double *ai, double *sc) {
  for (size_t q = 0; q < sz(p, p); ++q) sc[q] = a[q];
  if (!pev::chol(sc, p)) return false;
  pev::cinv(sc, p, ai);
  return true;
}

// The eigenvalues (ascending, into ev) and eigenvectors (columns of vec)
// of the symmetric p x p matrix a, by cyclic Jacobi rotations; a is
// overwritten. Returns false when 100 sweeps leave it off diagonal.
bool jacobi(double *a, int p, double *ev, double *vec) {
  const size_t sp = static_cast<size_t>(p);
  for (size_t q = 0; q < sp * sp; ++q) vec[q] = 0.0;
  for (size_t j = 0; j < sp; ++j) vec[j * sp + j] = 1.0;
  bool ok = false;
  for (int sw = 0; sw < 100 && !ok; ++sw) {
    double off = 0.0, tot = 0.0;
    for (size_t j = 0; j < sp; ++j) {
      for (size_t i = 0; i < sp; ++i) {
        tot += a[j * sp + i] * a[j * sp + i];
        if (i != j) off += a[j * sp + i] * a[j * sp + i];
      }
    }
    if (!(off > 1e-30 * tot)) {
      ok = true;
      break;
    }
    for (size_t r = 0; r + 1 < sp; ++r) {
      for (size_t c = r + 1; c < sp; ++c) {
        const double arc = a[c * sp + r];
        if (arc == 0.0) continue;
        const double th = (a[c * sp + c] - a[r * sp + r]) / (2 * arc);
        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double cs = 1 / std::sqrt(t * t + 1), sn = t * cs;
        for (size_t k = 0; k < sp; ++k) {
          const double akr = a[r * sp + k], akc = a[c * sp + k];
          a[r * sp + k] = cs * akr - sn * akc;
          a[c * sp + k] = sn * akr + cs * akc;
        }
        for (size_t k = 0; k < sp; ++k) {
          const double ark = a[k * sp + r], ack = a[k * sp + c];
          a[k * sp + r] = cs * ark - sn * ack;
          a[k * sp + c] = sn * ark + cs * ack;
        }
        for (size_t k = 0; k < sp; ++k) {
          const double vr = vec[r * sp + k], vc = vec[c * sp + k];
          vec[r * sp + k] = cs * vr - sn * vc;
          vec[c * sp + k] = sn * vr + cs * vc;
        }
      }
    }
  }
  // ascending order
  std::vector<size_t> ix(sp);
  for (size_t j = 0; j < sp; ++j) ix[j] = j;
  for (size_t j = 1; j < sp; ++j) {
    for (size_t i = j; i > 0 && a[ix[i] * sp + ix[i]] < a[ix[i - 1] * sp + ix[i - 1]]; --i) {
      const size_t t = ix[i];
      ix[i] = ix[i - 1];
      ix[i - 1] = t;
    }
  }
  std::vector<double> vv(vec, vec + sp * sp);
  for (size_t j = 0; j < sp; ++j) {
    ev[j] = a[ix[j] * sp + ix[j]];
    for (size_t k = 0; k < sp; ++k) vec[j * sp + k] = vv[ix[j] * sp + k];
  }
  return ok;
}

// The symmetric power s^e (e = 1/2 or -1/2) of the symmetric positive
// definite p x p matrix s into o; false when it is not positive definite.
bool spow(const double *s, int p, double e, double *o) {
  const size_t sp = static_cast<size_t>(p);
  std::vector<double> a(s, s + sp * sp), ev(sp), vec(sp * sp);
  if (!jacobi(a.data(), p, ev.data(), vec.data()) || !(ev[0] > 1e-12 * ev[sp - 1])) return false;
  for (size_t j = 0; j < sp; ++j) {
    for (size_t i = 0; i < sp; ++i) {
      double t = 0.0;
      for (size_t k = 0; k < sp; ++k) t += vec[k * sp + i] * std::pow(ev[k], e) * vec[k * sp + j];
      o[j * sp + i] = t;
    }
  }
  return true;
}

// The rk Wald statistic of Kleibergen and Paap (2006, equation 22, p.
// 104) that the first stages of pe > 1 endogenous columns have rank pe - 1
// (Corollary 1: chi-squared on q - pe + 1 degrees of freedom). Pi (pe x q)
// holds the first stages' instrument coefficients (pz, column j the jth
// stage), V its cluster-robust covariance from the residuals fr and the
// instruments rz less the controls, Q = rz'W rz; Theta = G Pi F' with F F'
// = Q and G'G = the inverse of the cross products of the endogenous
// columns less the controls (their normalisation, p. 103, with the sample
// moments in place of T^-1 sums, which only scales). Their Lambda_q (11)
// is the smallest singular value, pre- and post-multiplied by
// orthonormal matrices, so lambda = V2' Theta' u with u the left singular
// vector of the smallest singular value and V2 an orthonormal basis of
// the right singular vectors from the pe-th on; rk = lambda' Omega^-1
// lambda with Omega = (V2'F (x) u'G) V (V2'F (x) u'G)'. The statistic does
// not depend on which orthonormal bases are taken. NaN when a matrix is
// singular.
double rkwald(const double *fr, const double *pz, const double *rz, const double *qq,
              const double *ee, const double *ex, int n, int pe, int px, int q, const double *w,
              const int *cl, int ncl, const int *cl2, int ncl2, const int *c12, int n12, double *b,
              double *xi, double *sc) {
  const size_t sn = static_cast<size_t>(n), k = static_cast<size_t>(pe), m = static_cast<size_t>(q);
  const size_t km = k * m, d = m - k + 1;
  // the cross products of the endogenous columns less the controls
  std::vector<double> ey(sz(n, pe)), sy(k * k), qi(m * m), qs(m * m);
  for (int j = 0; j < pe; ++j) {
    const double *ej = ee + sz(n, j);
    if (px > 0) {
      if (pev::ols(ex, ej, n, px, w, b, xi, ey.data() + sz(n, j), sc) != CLAB_OK) return kutil::nan();
    } else {
      for (size_t i = 0; i < sn; ++i) ey[sz(n, j) + i] = ej[i];
    }
  }
  for (size_t j = 0; j < k; ++j) {
    for (size_t l = 0; l < k; ++l) {
      double t = 0.0;
      for (size_t i = 0; i < sn; ++i) t += (w ? w[i] : 1.0) * ey[j * sn + i] * ey[l * sn + i];
      sy[l * k + j] = t;
    }
  }
  if (!spdinv(qq, q, qi.data(), qs.data())) return kutil::nan();
  // the scores of vec(Pi), element (j, a) at a pe + j: w_i e_ij (Q^-1 rz_i)_a
  std::vector<double> xs(sn * km), one(sn, 1.0), vh(km * km), m2(km * km);
  for (size_t a = 0; a < m; ++a) {
    for (size_t i = 0; i < sn; ++i) {
      double t = 0.0;
      for (size_t c = 0; c < m; ++c) t += qi[c * m + a] * rz[c * sn + i];
      for (size_t j = 0; j < k; ++j) xs[(a * k + j) * sn + i] = t * fr[j * sn + i];
    }
  }
  const int g = cl2 && ncl2 < ncl ? ncl2 : ncl;
  const int big = ncl > ncl2 ? (ncl > n12 ? ncl : n12) : (ncl2 > n12 ? ncl2 : n12);
  std::vector<double> mk(static_cast<size_t>(big) * km);
  pev::meat(xs.data(), one.data(), w, n, static_cast<int>(km), cl, ncl, vh.data(), mk.data());
  if (cl2) {
    pev::meat(xs.data(), one.data(), w, n, static_cast<int>(km), cl2, ncl2, m2.data(), mk.data());
    for (size_t a = 0; a < km * km; ++a) vh[a] += m2[a];
    pev::meat(xs.data(), one.data(), w, n, static_cast<int>(km), c12, n12, m2.data(), mk.data());
    for (size_t a = 0; a < km * km; ++a) vh[a] -= m2[a];
  }
  for (size_t a = 0; a < km * km; ++a) vh[a] *= static_cast<double>(g) / (g - 1);
  // F = Q^(1/2), G = Sy^(-1/2), Theta = G Pi F
  std::vector<double> ff(m * m), gg(k * k), gp(k * m), th(k * m);
  if (!spow(qq, q, 0.5, ff.data()) || !spow(sy.data(), pe, -0.5, gg.data())) return kutil::nan();
  for (size_t j = 0; j < k; ++j) {
    for (size_t a = 0; a < m; ++a) {
      double t = 0.0;
      for (size_t l = 0; l < k; ++l) t += gg[l * k + j] * pz[l * m + a];
      gp[a * k + j] = t;
    }
  }
  for (size_t j = 0; j < k; ++j) {
    for (size_t a = 0; a < m; ++a) {
      double t = 0.0;
      for (size_t c = 0; c < m; ++c) t += gp[c * k + j] * ff[a * m + c];
      th[a * k + j] = t;
    }
  }
  std::vector<double> tt(k * k), ek(k), uk(k * k), ss(m * m), em(m), vm(m * m);
  for (size_t j = 0; j < k; ++j) {
    for (size_t l = 0; l < k; ++l) {
      double t = 0.0;
      for (size_t a = 0; a < m; ++a) t += th[a * k + j] * th[a * k + l];
      tt[l * k + j] = t;
    }
  }
  for (size_t a = 0; a < m; ++a) {
    for (size_t c = 0; c < m; ++c) {
      double t = 0.0;
      for (size_t j = 0; j < k; ++j) t += th[a * k + j] * th[c * k + j];
      ss[c * m + a] = t;
    }
  }
  if (!jacobi(tt.data(), pe, ek.data(), uk.data()) || !jacobi(ss.data(), q, em.data(), vm.data())) {
    return kutil::nan();
  }
  // u: the first (smallest) eigenvector of Theta Theta'; V2: the d
  // smallest of Theta'Theta
  const double *u = uk.data();
  std::vector<double> lam(d), gu(k), vf(d * m), mm(d * km), om(d * d), oi(d * d), os(d * d);
  for (size_t r = 0; r < d; ++r) {
    double t = 0.0;
    for (size_t a = 0; a < m; ++a) {
      double tu = 0.0;
      for (size_t j = 0; j < k; ++j) tu += th[a * k + j] * u[j];
      t += vm[r * m + a] * tu;
    }
    lam[r] = t;
    for (size_t a = 0; a < m; ++a) {
      double f = 0.0;
      for (size_t c = 0; c < m; ++c) f += vm[r * m + c] * ff[a * m + c];
      vf[a * d + r] = f;
    }
  }
  for (size_t j = 0; j < k; ++j) {
    double t = 0.0;
    for (size_t l = 0; l < k; ++l) t += u[l] * gg[j * k + l];
    gu[j] = t;
  }
  for (size_t r = 0; r < d; ++r) {
    for (size_t a = 0; a < m; ++a) {
      for (size_t j = 0; j < k; ++j) mm[(a * k + j) * d + r] = vf[a * d + r] * gu[j];
    }
  }
  for (size_t r = 0; r < d; ++r) {
    for (size_t c = 0; c < d; ++c) {
      double t = 0.0;
      for (size_t x = 0; x < km; ++x) {
        double s2 = 0.0;
        for (size_t y = 0; y < km; ++y) s2 += vh[y * km + x] * mm[y * d + c];
        t += mm[x * d + r] * s2;
      }
      om[c * d + r] = t;
    }
  }
  if (!spdinv(om.data(), static_cast<int>(d), oi.data(), os.data())) return kutil::nan();
  double rk = 0.0;
  for (size_t r = 0; r < d; ++r) {
    for (size_t c = 0; c < d; ++c) rk += lam[r] * oi[c * d + r] * lam[c];
  }
  return rk;
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
  return clab_ivfek(y, xe, pe, xx, px, zi, q, n, w, f1, n1, f2, n2, nullptr, nullptr, 0, cl, ncl,
                    cl2, ncl2, tol, beta, vc, fst, fres, res, sco, brd, st);
}

int clab_ivfek(const double *y, const double *xe, int pe, const double *xx, int px,
               const double *zi, int q, int n, const double *w, const int *f1, int n1,
               const int *f2, int n2, const int *fx, const int *nx, int mx, const int *cl,
               int ncl, const int *cl2, int ncl2, double tol, double *beta, double *vc,
               double *fst, double *fres, double *res, double *sco, double *brd, double *st) {
  const int p = pe + px, k = 1 + pe + px + q;
  if (n < 2 || pe < 0 || px < 0 || q < pe || p < 1 || n1 < 1 || (f2 && n2 < 1) ||
      !(tol > 0) || (cl && ncl < 2) || (cl2 && (!cl || ncl2 < 2)) ||
      !fecodes(n, f1, n1, f2, n2, fx, nx, mx)) {
    return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if ((cl && (cl[i] < 0 || cl[i] >= ncl)) || (cl2 && (cl2[i] < 0 || cl2[i] >= ncl2))) {
      return CLAB_EARG;
    }
    if (w && !(w[i] > 0 && std::isfinite(w[i]))) return CLAB_EDOM;
  }
  const size_t sn = static_cast<size_t>(n);
  const int m = q + px > p ? q + px : p;
  std::vector<double> work(ivwork(n, k, m, ngrp(n1, n2, nx, mx), 0));
  // v: y, xe, xx, zi, each demeaned
  double *v = work.data(), *r = v + sz(n, k), *u = r + sn;
  double *xi = u + sn, *xi2 = xi + sz(m, m), *b = xi2 + sz(m, m);
  double *sc = b + 2 * static_cast<size_t>(m);
  for (size_t i = 0; i < sn; ++i) v[i] = y[i];
  for (size_t a = 0; a < sz(n, pe); ++a) v[sn + a] = xe[a];
  for (size_t a = 0; a < sz(n, px); ++a) v[sn * (1 + pe) + a] = xx[a];
  for (size_t a = 0; a < sz(n, q); ++a) v[sn * (1 + pe + px) + a] = zi[a];
  const int sw = demean(v, n, k, f1, n1, f2, n2, w, tol, sc, fx, nx, mx);
  if (sw < 0) return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
  // a column the effects absorb keeps only rounding after the demeaning
  for (int j = 1; j < k; ++j) {
    const double *src = j < 1 + pe ? xe + sz(n, j - 1)
                        : j < 1 + pe + px ? xx + sz(n, j - 1 - pe) : zi + sz(n, j - 1 - pe - px);
    double a = 0.0, bv = 0.0;
    for (size_t i = 0; i < sn; ++i) {
      a += src[i] * src[i];
      bv += v[sz(n, j) + i] * v[sz(n, j) + i];
    }
    if (!(bv > 1e-16 * a)) return CLAB_EDOM;
  }
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
  // the instruments less their fit on the controls (all demeaned), and
  // their weighted cross products Q
  std::vector<double> rz, qq;
  if (cl) {
    rz.resize(sz(n, q));
    qq.resize(sz(q, q));
    for (int a = 0; a < q; ++a) {
      const double *za = fz.data() + sz(n, a);
      if (px > 0) {
        if (pev::ols(ex, za, n, px, w, b, xi, rz.data() + sz(n, a), sc) != CLAB_OK) return CLAB_EDOM;
      } else {
        for (size_t i = 0; i < sn; ++i) rz[sz(n, a) + i] = za[i];
      }
    }
    for (int a = 0; a < q; ++a) {
      for (int c = 0; c < q; ++c) {
        kutil::ksum t;
        for (size_t i = 0; i < sn; ++i) t.add((w ? w[i] : 1.0) * rz[sz(n, a) + i] * rz[sz(n, c) + i]);
        qq[sz(c, q) + static_cast<size_t>(a)] = t.val();
      }
    }
  }
  // each first stage's residuals and instrument coefficients, for the
  // rank statistic of several endogenous columns
  std::vector<double> fr(sz(n, pe)), pz(sz(q, pe));
  for (int j = 0; j < pe; ++j) {
    const double *ej = ee + sz(n, j);
    if (pev::ols(fz.data(), ej, n, q + px, w, b, xi, r, sc) != CLAB_OK) return CLAB_EDOM;
    for (size_t i = 0; i < sn; ++i) {
      h[sz(n, j) + i] = ej[i] - r[i];
      fr[sz(n, j) + i] = r[i];
    }
    for (int a = 0; a < q; ++a) pz[sz(j, q) + static_cast<size_t>(a)] = b[a];
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
      double eff = kutil::nan(), keff = kutil::nan();
      if (pe == 1) {
        double num = 0.0, den = 0.0;
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) {
            num += bz[a] * qq[sz(c, q) + static_cast<size_t>(a)] * bz[c];
            den += vz[sz(c, q) + static_cast<size_t>(a)] * qq[sz(a, q) + static_cast<size_t>(c)];
          }
        }
        if (den > 0) eff = num / den;
        // the effective degrees of freedom of the simplified test,
        // Montiel Olea and Pflueger (2013, equation 4, p. 360) with x =
        // 1 / tau = 10: tr(M)^2 (1 + 2x) / (tr(M M) + 2x tr(M) maxeval(M)),
        // M = Q V (similar to their normalised W2, so the same traces and
        // eigenvalues)
        std::vector<double> mq(sz(q, q)), ev(static_cast<size_t>(q)), vec(sz(q, q));
        for (int a = 0; a < q; ++a) {
          for (int c = 0; c < q; ++c) {
            double t = 0.0;
            for (int e = 0; e < q; ++e) {
              t += qq[sz(e, q) + static_cast<size_t>(a)] * vz[sz(c, q) + static_cast<size_t>(e)];
            }
            mq[sz(c, q) + static_cast<size_t>(a)] = t;
          }
        }
        double tr = 0.0, tr2 = 0.0;
        for (int a = 0; a < q; ++a) {
          tr += mq[sz(a, q) + static_cast<size_t>(a)];
          for (int c = 0; c < q; ++c) {
            tr2 += mq[sz(c, q) + static_cast<size_t>(a)] * mq[sz(a, q) + static_cast<size_t>(c)];
          }
        }
        // the largest eigenvalue of Q V: that of L'VL, Q = LL'
        std::vector<double> lq(qq), lvl(sz(q, q));
        if (pev::chol(lq.data(), q)) {
          for (int a = 0; a < q; ++a) {
            for (int c = 0; c < q; ++c) {
              double t = 0.0;
              for (int e = a; e < q; ++e) {
                for (int f = c; f < q; ++f) {
                  t += lq[sz(a, q) + static_cast<size_t>(e)] * vz[sz(f, q) + static_cast<size_t>(e)] *
                       lq[sz(c, q) + static_cast<size_t>(f)];
                }
              }
              lvl[sz(c, q) + static_cast<size_t>(a)] = t;
            }
          }
          if (jacobi(lvl.data(), q, ev.data(), vec.data()) && tr > 0) {
            const double x = 10.0;
            keff = tr * tr * (1 + 2 * x) / (tr2 + 2 * x * tr * ev[static_cast<size_t>(q) - 1]);
          }
        }
      }
      fst[2 * pe + j] = eff;
      fst[3 * pe + j] = keff;
    }
  }
  if (cl && pe > 1) {
    const double rk = rkwald(fr.data(), pz.data(), rz.data(), qq.data(), ee, ex, n, pe, px, q, w, cl,
                             ncl, cl2, ncl2, c12.data(), n12, b, xi, sc);
    for (int j = 0; j < pe; ++j) fst[pe + j] = rk / q;
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
  double s0[7];
  int rc = one(all, nullptr, cl, ncl, beta, fst, s0);
  if (rc != CLAB_OK) return rc;
  st[0] = s0[0];
  st[1] = s0[1];
  int bad = 0;
  std::vector<double> fw(static_cast<size_t>(n)), sep(static_cast<size_t>(n1 + n2));
  std::vector<double> fsb(4 * static_cast<size_t>(pe)), sb(7), bb(static_cast<size_t>(p));
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

int clab_expgmm(const double *y, const double *x, int k, const double *z, int s, int n,
                const double *off, const double *b0, const int *cl, int ncl, int err, int step,
                int maxit, double tol, double *beta, double *vc, double *st) {
  if (n < 2 || k < 1 || s < k || !cl || ncl < 2 || err < 0 || err > 1 || step < 1 || step > 3 ||
      maxit < 1 || !(tol > 0)) {
    return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if (cl[i] < 0 || cl[i] >= ncl) return CLAB_EARG;
    if (!(y[i] >= 0) || !std::isfinite(y[i]) || (off && !std::isfinite(off[i]))) return CLAB_EDOM;
  }
  for (size_t a = 0; a < sz(n, k); ++a) {
    if (!std::isfinite(x[a])) return CLAB_EDOM;
  }
  for (size_t a = 0; a < sz(n, s); ++a) {
    if (!std::isfinite(z[a])) return CLAB_EDOM;
  }
  const size_t sn = static_cast<size_t>(n), sk = static_cast<size_t>(k);
  const size_t ss = static_cast<size_t>(s);
  std::vector<double> u(sn), d(sz(n, k)), g(ss), gg(sz(s, k)), w(sz(s, s)), sm(sz(s, s));
  std::vector<double> wg(sz(s, k)), a(sz(k, k)), ai(sz(k, k)), h(sk), dl(sk), bt(b0, b0 + k);
  std::vector<double> tb(sk), bp(sk), sc(ss > sk ? sz(s, s) : sz(k, k));
  std::vector<double> mk(sz(ncl, s)), mm(sz(k, k)), wsw(sz(s, k));
  // u (n) and its derivative d (n x k) at b; false when not finite. mult:
  // u = y exp(-x'b) - 1 (Mullahy 1997, the residual of IV-T); add: u = y -
  // exp(x'b) (Windmeijer and Santos Silva 1997, equation 2)
  auto res = [&](const double *b) -> bool {
    for (size_t i = 0; i < sn; ++i) {
      double eta = off ? off[i] : 0.0;
      for (size_t j = 0; j < sk; ++j) eta += x[j * sn + i] * b[j];
      double ui, di;
      if (err == 0) {
        const double e = std::exp(-eta);
        ui = y[i] * e - 1.0;
        di = -y[i] * e;
      } else {
        const double m = std::exp(eta);
        ui = y[i] - m;
        di = -m;
      }
      if (!std::isfinite(ui) || !std::isfinite(di)) return false;
      u[i] = ui;
      for (size_t j = 0; j < sk; ++j) d[j * sn + i] = di * x[j * sn + i];
    }
    return true;
  };
  // g = Z'u (s) and G = Z'D (s x k), sums over the rows
  auto mom = [&]() {
    for (size_t c = 0; c < ss; ++c) {
      const double *zc = z + c * sn;
      double t = 0.0;
      for (size_t i = 0; i < sn; ++i) t += zc[i] * u[i];
      g[c] = t;
      for (size_t j = 0; j < sk; ++j) {
        double v = 0.0;
        for (size_t i = 0; i < sn; ++i) v += zc[i] * d[j * sn + i];
        gg[j * ss + c] = v;
      }
    }
  };
  auto obj = [&]() {
    double q = 0.0;
    for (size_t c = 0; c < ss; ++c) {
      for (size_t e = 0; e < ss; ++e) q += g[e] * w[c * ss + e] * g[c];
    }
    return q;
  };
  // A = G'WG into a, its inverse into ai, h = G'Wg; false when singular
  auto norm = [&]() -> bool {
    for (size_t j = 0; j < sk; ++j) {
      for (size_t c = 0; c < ss; ++c) {
        double v = 0.0;
        for (size_t e = 0; e < ss; ++e) v += w[e * ss + c] * gg[j * ss + e];
        wg[j * ss + c] = v;
      }
    }
    for (size_t j = 0; j < sk; ++j) {
      double t = 0.0;
      for (size_t c = 0; c < ss; ++c) t += wg[j * ss + c] * g[c];
      h[j] = t;
      for (size_t l = 0; l < sk; ++l) {
        double v = 0.0;
        for (size_t c = 0; c < ss; ++c) v += gg[j * ss + c] * wg[l * ss + c];
        a[l * sk + j] = v;
      }
    }
    return spdinv(a.data(), k, ai.data(), sc.data());
  };
  // the inverse of the cluster sum of Z_g'u_g (Z_g'u_g)' at the current u
  // into w; false when singular (fewer clusters than moments, say)
  auto wopt = [&]() -> bool {
    pev::meat(z, u.data(), nullptr, n, s, cl, ncl, sm.data(), mk.data());
    return spdinv(sm.data(), s, w.data(), sc.data());
  };
  // Gauss-Newton on g'Wg from bt, halving a step that does not lower it;
  // returns the iterations, -1 when the start is not finite, -2 when
  // G'WG is singular, -3 on an interrupt
  int conv = 0;
  auto fit = [&]() -> int {
    conv = 0;
    if (!res(bt.data())) return -1;
    mom();
    double q = obj();
    for (int it = 1; it <= maxit; ++it) {
      if (kutil::stop(static_cast<double>(n) * s * (k + 1))) return -3;
      if (!norm()) return -2;
      for (size_t j = 0; j < sk; ++j) {
        double v = 0.0;
        for (size_t l = 0; l < sk; ++l) v += ai[l * sk + j] * h[l];
        dl[j] = -v;
      }
      double lam = 1.0;
      bool ok = false;
      for (int hv = 0; hv <= 50 && !ok; ++hv, lam /= 2) {
        for (size_t j = 0; j < sk; ++j) tb[j] = bt[j] + lam * dl[j];
        if (!res(tb.data())) continue;
        mom();
        const double qn = obj();
        if (qn <= q) {
          ok = true;
          q = qn;
        }
      }
      if (!ok) {
        // no step lowers the criterion: a minimum to working precision
        res(bt.data());
        mom();
        conv = 1;
        return it;
      }
      double ch = 0.0;
      for (size_t j = 0; j < sk; ++j) {
        ch = std::fmax(ch, std::fabs(tb[j] - bt[j]) / (1.0 + std::fabs(bt[j])));
        bt[j] = tb[j];
      }
      if (ch < tol) {
        conv = 1;
        return it;
      }
    }
    return maxit;
  };
  auto fail = [](int r) { return r == -3 ? CLAB_EINTR : CLAB_EDOM; };
  // the first step: the nonlinear instrumental-variables weight (Z'Z)^-1
  for (size_t c = 0; c < ss; ++c) {
    for (size_t e = 0; e < ss; ++e) {
      double v = 0.0;
      for (size_t i = 0; i < sn; ++i) v += z[c * sn + i] * z[e * sn + i];
      sm[c * ss + e] = v;
    }
  }
  if (!spdinv(sm.data(), s, w.data(), sc.data())) return CLAB_EDOM;
  int it = fit();
  if (it < 0) return fail(it);
  int nw = 1;
  while (step >= 2) {
    // the optimal weight from the last estimate's residuals
    if (!wopt()) return CLAB_EDOM;
    bp = bt;
    it = fit();
    if (it < 0) return fail(it);
    ++nw;
    if (step == 2) break;
    double ch = 0.0;
    for (size_t j = 0; j < sk; ++j) {
      ch = std::fmax(ch, std::fabs(bt[j] - bp[j]) / (1.0 + std::fabs(bp[j])));
    }
    if (ch < tol) break;
    if (nw > maxit) {
      conv = 0;
      break;
    }
  }
  if (!norm()) return CLAB_EDOM;
  for (size_t q = 0; q < sz(k, k); ++q) vc[q] = 0.0;
  if (step == 1) {
    // the sandwich (G'WG)^-1 G'W S W G (G'WG)^-1, S at the estimate
    pev::meat(z, u.data(), nullptr, n, s, cl, ncl, sm.data(), mk.data());
    for (size_t j = 0; j < sk; ++j) {
      for (size_t c = 0; c < ss; ++c) {
        double v = 0.0;
        for (size_t e = 0; e < ss; ++e) v += sm[e * ss + c] * wg[j * ss + e];
        wsw[j * ss + c] = v;
      }
    }
    for (size_t j = 0; j < sk; ++j) {
      for (size_t l = 0; l < sk; ++l) {
        double v = 0.0;
        for (size_t c = 0; c < ss; ++c) v += wg[j * ss + c] * wsw[l * ss + c];
        mm[l * sk + j] = v;
      }
    }
    pev::sand(ai.data(), mm.data(), k, 1.0, vc, sc.data());
  } else {
    // the efficient (G'WG)^-1 with W the inverse of S at the last step
    for (size_t q = 0; q < sz(k, k); ++q) vc[q] = ai[q];
  }
  for (size_t j = 0; j < sk; ++j) beta[j] = bt[j];
  // Hansen's J: the criterion at the efficient estimate (Mullahy 1997,
  // p. 590), chi-squared on s - k degrees of freedom
  st[0] = step >= 2 && s > k ? obj() : kutil::nan();
  st[1] = s - k;
  st[2] = it;
  st[3] = conv;
  st[4] = nw;
  st[5] = obj();
  return CLAB_OK;
}

int clab_ivoid(const double *y, const double *xe, int pe, const double *xx, int px,
               const double *zi, int q, int n, const int *f1, int n1, const int *f2, int n2,
               const int *fx, const int *nx, int mx, const int *cl, int ncl, const int *cl2,
               int ncl2, double tol, double *beta, double *st) {
  const int p = pe + px, m = q + px, k = 1 + pe + px + q;
  if (n < 2 || pe < 1 || px < 0 || q < pe || !(tol > 0) || !cl || ncl < 2 ||
      (cl2 && ncl2 < 2) || !fecodes(n, f1, n1, f2, n2, fx, nx, mx)) {
    return CLAB_EARG;
  }
  for (int i = 0; i < n; ++i) {
    if (cl[i] < 0 || cl[i] >= ncl || (cl2 && (cl2[i] < 0 || cl2[i] >= ncl2))) return CLAB_EARG;
  }
  const size_t sn = static_cast<size_t>(n), sp = static_cast<size_t>(p);
  const size_t smm = static_cast<size_t>(m);
  // v: y, xe, xx, zi, each demeaned
  std::vector<double> v(sz(n, k));
  for (size_t i = 0; i < sn; ++i) v[i] = y[i];
  for (size_t a = 0; a < sz(n, pe); ++a) v[sn + a] = xe[a];
  for (size_t a = 0; a < sz(n, px); ++a) v[sn * (1 + pe) + a] = xx[a];
  for (size_t a = 0; a < sz(n, q); ++a) v[sn * (1 + pe + px) + a] = zi[a];
  std::vector<double> wk(4 * static_cast<size_t>(ngrp(n1, n2, nx, mx)));
  if (demean(v.data(), n, k, f1, n1, f2, n2, nullptr, tol, wk.data(), fx, nx, mx) < 0) {
    return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
  }
  // X = (xe, xx) and Z = (zi, xx), demeaned
  std::vector<double> xm(sz(n, p)), zm(sz(n, m));
  for (size_t a = 0; a < sz(n, p); ++a) xm[a] = v[sn + a];
  for (size_t a = 0; a < sz(n, q); ++a) zm[a] = v[sn * (1 + pe + px) + a];
  for (size_t a = 0; a < sz(n, px); ++a) zm[sz(n, q) + a] = v[sn * (1 + pe) + a];
  const double *ye = v.data();
  std::vector<double> zz(sz(m, m)), zzi(sz(m, m)), zx(sz(m, p)), zy(smm), zu(smm);
  std::vector<double> sc(sz(m > p ? m : p, m > p ? m : p)), wz(sz(m, p)), aa(sz(p, p));
  std::vector<double> ai(sz(p, p)), hh(sp), u(sn), wm(sz(m, m)), sm(sz(m, m)), sm2(sz(m, m));
  for (size_t c = 0; c < smm; ++c) {
    for (size_t e = 0; e < smm; ++e) {
      double t = 0.0;
      for (size_t i = 0; i < sn; ++i) t += zm[c * sn + i] * zm[e * sn + i];
      zz[e * smm + c] = t;
    }
    double t = 0.0;
    for (size_t i = 0; i < sn; ++i) t += zm[c * sn + i] * ye[i];
    zy[c] = t;
    for (size_t j = 0; j < sp; ++j) {
      double r = 0.0;
      for (size_t i = 0; i < sn; ++i) r += zm[c * sn + i] * xm[j * sn + i];
      zx[j * smm + c] = r;
    }
  }
  if (!spdinv(zz.data(), m, zzi.data(), sc.data())) return CLAB_EDOM;
  // the linear GMM estimate with weight wt (m x m) into b; residuals into
  // u and Z'u into zu; false when X'Z W Z'X is singular
  auto lgmm = [&](const double *wt, double *b) -> bool {
    for (size_t j = 0; j < sp; ++j) {
      for (size_t c = 0; c < smm; ++c) {
        double t = 0.0;
        for (size_t e = 0; e < smm; ++e) t += wt[e * smm + c] * zx[j * smm + e];
        wz[j * smm + c] = t;
      }
    }
    for (size_t j = 0; j < sp; ++j) {
      double t = 0.0;
      for (size_t c = 0; c < smm; ++c) t += wz[j * smm + c] * zy[c];
      hh[j] = t;
      for (size_t l = 0; l < sp; ++l) {
        double r = 0.0;
        for (size_t c = 0; c < smm; ++c) r += zx[j * smm + c] * wz[l * smm + c];
        aa[l * sp + j] = r;
      }
    }
    if (!spdinv(aa.data(), p, ai.data(), sc.data())) return false;
    for (size_t j = 0; j < sp; ++j) {
      double t = 0.0;
      for (size_t l = 0; l < sp; ++l) t += ai[l * sp + j] * hh[l];
      b[j] = t;
    }
    for (size_t i = 0; i < sn; ++i) {
      double f = ye[i];
      for (size_t j = 0; j < sp; ++j) f -= xm[j * sn + i] * b[j];
      u[i] = f;
    }
    for (size_t c = 0; c < smm; ++c) {
      double t = 0.0;
      for (size_t i = 0; i < sn; ++i) t += zm[c * sn + i] * u[i];
      zu[c] = t;
    }
    return true;
  };
  auto quad = [&](const double *wt) {
    double t = 0.0;
    for (size_t c = 0; c < smm; ++c) {
      for (size_t e = 0; e < smm; ++e) t += zu[e] * wt[c * smm + e] * zu[c];
    }
    return t;
  };
  // two-stage least squares (weight (Z'Z)^-1) and the Sargan statistic
  // n u'Pz u / u'u
  std::vector<double> b1(sp);
  if (!lgmm(zzi.data(), b1.data())) return CLAB_EDOM;
  double uu = 0.0;
  for (size_t i = 0; i < sn; ++i) uu += u[i] * u[i];
  st[0] = uu > 0 ? n * quad(zzi.data()) / uu : kutil::nan();
  // the cluster-robust S of Z'u at the 2SLS residuals, one-way or
  // V1 + V2 - V12, then the two-step estimate and Hansen's J
  std::vector<double> mk(sz(ncl > ncl2 ? ncl : ncl2, m));
  pev::meat(zm.data(), u.data(), nullptr, n, m, cl, ncl, sm.data(), mk.data());
  if (cl2) {
    std::vector<int> c12(sn), pw(sz(ncl, ncl2));
    const int n12 = pev::pairs(cl, cl2, n, ncl, ncl2, c12.data(), pw.data());
    std::vector<double> m12(sz(n12, m));
    pev::meat(zm.data(), u.data(), nullptr, n, m, cl2, ncl2, sm2.data(), mk.data());
    for (size_t a = 0; a < sz(m, m); ++a) sm[a] += sm2[a];
    pev::meat(zm.data(), u.data(), nullptr, n, m, c12.data(), n12, sm2.data(), m12.data());
    for (size_t a = 0; a < sz(m, m); ++a) sm[a] -= sm2[a];
  }
  st[1] = kutil::nan();
  st[2] = q - pe;
  for (size_t j = 0; j < sp; ++j) beta[j] = kutil::nan();
  if (spdinv(sm.data(), m, wm.data(), sc.data()) && lgmm(wm.data(), beta)) {
    st[1] = q > pe ? quad(wm.data()) : kutil::nan();
  }
  return CLAB_OK;
}

int clab_demean(double *v, int n, int k, const int *f1, int n1, const int *f2, int n2,
                const int *fx, const int *nx, int mx, double tol, double *work) {
  if (n < 1 || k < 0 || !(tol > 0) || !fecodes(n, f1, n1, f2, n2, fx, nx, mx)) return CLAB_EARG;
  if (demean(v, n, k, f1, n1, f2, n2, nullptr, tol, work, fx, nx, mx) < 0) {
    return kutil::stop(0) ? CLAB_EINTR : CLAB_EDOM;
  }
  return CLAB_OK;
}

int clab_windbin(const double *dir, int n, double width, double start, const int *area, int na,
                 int base, int edge, double *out) {
  if (n < 1 || !(width > 0) || !(width < 360) || !std::isfinite(start) || na < 1 || edge < 0 ||
      edge > 1) {
    return CLAB_EARG;
  }
  const double nbd = 360 / width;
  // the sectors, and the columns (areas times sectors) counted in int
  if (nbd > kutil::maxcell || static_cast<double>(na) * (std::nearbyint(nbd) - 1) > 2147483647.0) return CLAB_EARG;
  const int nb = static_cast<int>(std::nearbyint(nbd));
  if (std::fabs(nbd - nb) > 1e-9 || nb < 2 || base < 0 || base >= nb) return CLAB_EARG;
  const int k = nb - 1;
  for (size_t a = 0; a < sz(n, na * k); ++a) out[a] = 0.0;
  for (int i = 0; i < n; ++i) {
    const int ar = area ? area[i] : 0;
    if (ar < 0 || ar >= na) return CLAB_EARG;
    if (!std::isfinite(dir[i])) {
      for (std::size_t c = 0; c < static_cast<std::size_t>(na) * static_cast<std::size_t>(k); ++c) out[sz(c, n) + static_cast<size_t>(i)] = kutil::nan();
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

int clab_runlen(const int *u, const double *t, const double *e, int n, double step, int *run,
                int *ord) {
  if (n < 1 || !(step > 0)) return CLAB_EARG;
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(t[i])) return CLAB_EDOM;
    ord[i] = i;
  }
  auto less = [&](int a, int b) { return u[a] < u[b] || (u[a] == u[b] && t[a] < t[b]); };
  kutil::hsort(ord, n, less);
  for (int q = 0; q < n; ++q) {
    const int i = ord[q];
    if (q > 0 && u[ord[q - 1]] == u[i] && t[ord[q - 1]] == t[i]) return CLAB_EDOM;
    if (!std::isfinite(e[i])) {
      run[i] = -1;
      continue;
    }
    if (!(e[i] > 0)) {
      run[i] = 0;
      continue;
    }
    // the run goes on from the row of the same unit one step earlier
    const int j = q > 0 ? ord[q - 1] : -1;
    const bool on = j >= 0 && u[j] == u[i] && t[i] - t[j] == step && run[j] > 0;
    run[i] = on ? run[j] + 1 : 1;
  }
  return CLAB_OK;
}

}  // extern "C"
