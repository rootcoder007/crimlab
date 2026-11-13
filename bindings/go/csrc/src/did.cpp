// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Difference in differences on unit-by-period panels: the within
// transformation and weighted least squares with cluster-robust variances
// (pev.h), the two-way fixed-effects regression (clab_didfe), the
// Goodman-Bacon decomposition of its estimate (pev::bacon), and the
// nearest-neighbour matching of treated units on their outcomes before
// treatment (pev::mtch).
#include <cmath>
#include <cstddef>

#include "crimlab.h"
#include "kutil.h"
#include "pev.h"

namespace pev {

namespace {

inline std::size_t z(int i) { return static_cast<std::size_t>(i); }

// weighted means of c by the codes g (ng groups, weight totals sw) into m
void gmean(const double *c, const double *w, const int *g, int n, int ng, const double *sw,
           double *m) {
  for (int a = 0; a < ng; ++a) m[a] = 0.0;
  for (int i = 0; i < n; ++i) m[g[i]] += (w ? w[i] : 1.0) * c[i];
  for (int a = 0; a < ng; ++a) m[a] /= sw[a];
}

}  // namespace

int within(double *v, int n, int k, const int *u, int nu, const int *t, int nt, int two,
           const double *w, double tol, int maxit, double *work) {
  double *su = work, *st = work + nu, *mu = work + nu + nt, *mt = work + 2 * nu + nt;
  for (int a = 0; a < nu; ++a) su[a] = 0.0;
  for (int a = 0; a < nt; ++a) st[a] = 0.0;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[i] : 1.0;
    if (two) su[u[i]] += wi;
    st[t[i]] += wi;
  }
  int most = two ? 0 : 1;
  for (int j = 0; j < k; ++j) {
    double *c = v + z(j) * z(n);
    if (!two) {
      gmean(c, w, t, n, nt, st, mt);
      for (int i = 0; i < n; ++i) c[i] -= mt[t[i]];
      continue;
    }
    double big = 0.0;
    for (int i = 0; i < n; ++i) big = std::fmax(big, std::fabs(c[i]));
    const double lim = tol * (big > 0.0 ? big : 1.0);
    int s = 0;
    for (;;) {
      gmean(c, w, u, n, nu, su, mu);
      double left = 0.0;
      for (int a = 0; a < nu; ++a) left = std::fmax(left, std::fabs(mu[a]));
      if (s > 0 && left <= lim) break;
      if (s == maxit) return -1;
      for (int i = 0; i < n; ++i) c[i] -= mu[u[i]];
      gmean(c, w, t, n, nt, st, mt);
      for (int i = 0; i < n; ++i) c[i] -= mt[t[i]];
      ++s;
    }
    if (s > most) most = s;
  }
  return most;
}

bool chol(double *a, int p) {
  double big = 0.0;
  for (int j = 0; j < p; ++j) big = std::fmax(big, a[z(j) * z(p) + z(j)]);
  for (int j = 0; j < p; ++j) {
    double d = a[z(j) * z(p) + z(j)];
    for (int k = 0; k < j; ++k) d -= a[z(k) * z(p) + z(j)] * a[z(k) * z(p) + z(j)];
    if (!(d > 1e-12 * big)) return false;
    d = std::sqrt(d);
    a[z(j) * z(p) + z(j)] = d;
    for (int i = j + 1; i < p; ++i) {
      double s = a[z(j) * z(p) + z(i)];
      for (int k = 0; k < j; ++k) s -= a[z(k) * z(p) + z(i)] * a[z(k) * z(p) + z(j)];
      a[z(j) * z(p) + z(i)] = s / d;
    }
  }
  return true;
}

void cinv(const double *l, int p, double *ai) {
  // column by column: solve l l' x = e_c
  for (int c = 0; c < p; ++c) {
    double *x = ai + z(c) * z(p);
    for (int i = 0; i < p; ++i) {
      double s = i == c ? 1.0 : 0.0;
      for (int k = 0; k < i; ++k) s -= l[z(k) * z(p) + z(i)] * x[k];
      x[i] = s / l[z(i) * z(p) + z(i)];
    }
    for (int i = p - 1; i >= 0; --i) {
      double s = x[i];
      for (int k = i + 1; k < p; ++k) s -= l[z(i) * z(p) + z(k)] * x[k];
      x[i] = s / l[z(i) * z(p) + z(i)];
    }
  }
}

int ols(const double *x, const double *y, int n, int p, const double *w, double *beta,
        double *xi, double *r, double *work) {
  double *a = work, *b = work + z(p) * z(p);
  for (int j = 0; j < p; ++j) {
    const double *xj = x + z(j) * z(n);
    for (int k = j; k < p; ++k) {
      const double *xk = x + z(k) * z(n);
      kutil::ksum s;
      for (int i = 0; i < n; ++i) s.add((w ? w[i] : 1.0) * xj[i] * xk[i]);
      a[z(j) * z(p) + z(k)] = a[z(k) * z(p) + z(j)] = s.val();
    }
    kutil::ksum s;
    for (int i = 0; i < n; ++i) s.add((w ? w[i] : 1.0) * xj[i] * y[i]);
    b[j] = s.val();
  }
  if (!chol(a, p)) return CLAB_EDOM;
  cinv(a, p, xi);
  for (int j = 0; j < p; ++j) {
    double s = 0.0;
    for (int k = 0; k < p; ++k) s += xi[z(k) * z(p) + z(j)] * b[k];
    beta[j] = s;
  }
  for (int i = 0; i < n; ++i) {
    double f = 0.0;
    for (int j = 0; j < p; ++j) f += x[z(j) * z(n) + z(i)] * beta[j];
    r[i] = y[i] - f;
  }
  return CLAB_OK;
}

void meat(const double *x, const double *r, const double *w, int n, int p, const int *cl,
          int ncl, double *m, double *work) {
  for (std::size_t q = 0; q < z(ncl) * z(p); ++q) work[q] = 0.0;
  for (int j = 0; j < p; ++j) {
    for (int i = 0; i < n; ++i) {
      work[z(cl[i]) * z(p) + z(j)] += (w ? w[i] : 1.0) * r[i] * x[z(j) * z(n) + z(i)];
    }
  }
  for (std::size_t q = 0; q < z(p) * z(p); ++q) m[q] = 0.0;
  for (int g = 0; g < ncl; ++g) {
    const double *s = work + z(g) * z(p);
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < p; ++k) m[z(k) * z(p) + z(j)] += s[j] * s[k];
    }
  }
}

void sand(const double *xi, const double *m, int p, double c, double *v, double *work) {
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      double s = 0.0;
      for (int a = 0; a < p; ++a) s += xi[z(a) * z(p) + z(j)] * m[z(k) * z(p) + z(a)];
      work[z(k) * z(p) + z(j)] = s;
    }
  }
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < p; ++k) {
      double s = 0.0;
      for (int a = 0; a < p; ++a) s += work[z(a) * z(p) + z(j)] * xi[z(k) * z(p) + z(a)];
      v[z(k) * z(p) + z(j)] += c * s;
    }
  }
}

int crv(const double *x, const double *r, const double *w, int n, int p, const double *xi,
        const int *cl, int ncl, const int *cl2, int ncl2, const int *c12, int n12, double kk,
        double *v, double *work) {
  const int g = cl2 && ncl2 < ncl ? ncl2 : ncl;
  const double c = static_cast<double>(g) / (g - 1) * (n - 1.0) / (n - kk);
  double *m = work, *sw = work + z(p) * z(p), *sc = sw + z(p) * z(p);
  for (std::size_t q = 0; q < z(p) * z(p); ++q) v[q] = 0.0;
  meat(x, r, w, n, p, cl, ncl, m, sc);
  sand(xi, m, p, c, v, sw);
  if (cl2) {
    meat(x, r, w, n, p, cl2, ncl2, m, sc);
    sand(xi, m, p, c, v, sw);
    meat(x, r, w, n, p, c12, n12, m, sc);
    sand(xi, m, p, -c, v, sw);
  }
  return g;
}

int pairs(const int *a, const int *b, int n, int na, int nb, int *c, int *iwork) {
  for (std::size_t q = 0; q < z(na) * z(nb); ++q) iwork[q] = -1;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    int &s = iwork[z(a[i]) * z(nb) + z(b[i])];
    if (s < 0) s = k++;
    c[i] = s;
  }
  return k;
}

int bacon(const double *yw, int nu, int nt, const int *g, int *ty, int *tr, int *co,
          double *est, double *wt, double *work, int *iwork) {
  // the timing groups: distinct g in 0..nt-1, then the never treated (g >= nt)
  int *gl = iwork, *gi = iwork + nu;  // group values; each unit's group
  int ng = 0;
  bool never = false;
  for (int i = 0; i < nu; ++i) {
    if (g[i] >= nt) {
      never = true;
      continue;
    }
    int a = 0;
    while (a < ng && gl[a] != g[i]) ++a;
    if (a == ng) gl[ng++] = g[i];
  }
  // sort the group values (few of them)
  for (int a = 1; a < ng; ++a) {
    for (int b = a; b > 0 && gl[b - 1] > gl[b]; --b) {
      const int s = gl[b];
      gl[b] = gl[b - 1];
      gl[b - 1] = s;
    }
  }
  const int nk = ng + (never ? 1 : 0);
  for (int i = 0; i < nu; ++i) {
    int a = ng;
    if (g[i] < nt) {
      a = 0;
      while (gl[a] != g[i]) ++a;
    }
    gi[i] = a;
  }
  // group sizes and mean outcome by period: m (nk x nt)
  double *nn = work, *m = work + nk;
  for (int a = 0; a < nk; ++a) nn[a] = 0.0;
  for (std::size_t q = 0; q < z(nk) * z(nt); ++q) m[q] = 0.0;
  for (int i = 0; i < nu; ++i) {
    nn[gi[i]] += 1.0;
    for (int s = 0; s < nt; ++s) m[z(gi[i]) * z(nt) + z(s)] += yw[z(s) * z(nu) + z(i)];
  }
  for (int a = 0; a < nk; ++a) {
    for (int s = 0; s < nt; ++s) m[z(a) * z(nt) + z(s)] /= nn[a];
  }
  // mean of group a over periods [lo, hi)
  auto avg = [&](int a, int lo, int hi) {
    double s = 0.0;
    for (int q = lo; q < hi; ++q) s += m[z(a) * z(nt) + z(q)];
    return s / (hi - lo);
  };
  // the 2 x 2 estimate: treated a, control b, pre [lo, mid), post [mid, hi)
  auto dd = [&](int a, int b, int lo, int mid, int hi) {
    return (avg(a, mid, hi) - avg(a, lo, mid)) - (avg(b, mid, hi) - avg(b, lo, mid));
  };
  auto share = [&](int gv) { return static_cast<double>(nt - gv) / nt; };
  int r = 0;
  if (never) {
    for (int a = 0; a < ng; ++a) {
      if (gl[a] == 0) continue;  // always treated: no variation against never
      const double d = share(gl[a]);
      ty[r] = 0;
      tr[r] = gl[a];
      co[r] = nt;
      est[r] = dd(a, ng, 0, gl[a], nt);
      wt[r] = nn[a] * nn[ng] * d * (1.0 - d);
      ++r;
    }
  }
  for (int a = 0; a < ng; ++a) {
    for (int b = a + 1; b < ng; ++b) {
      const double dk = share(gl[a]), dl = share(gl[b]);
      if (gl[a] > 0) {
        // earlier a against later b, before b starts
        ty[r] = 1;
        tr[r] = gl[a];
        co[r] = gl[b];
        est[r] = dd(a, b, 0, gl[a], gl[b]);
        wt[r] = nn[a] * nn[b] * (dk - dl) * (1.0 - dk);
        ++r;
      }
      // later b against earlier a, once a has started
      ty[r] = 2;
      tr[r] = gl[b];
      co[r] = gl[a];
      est[r] = dd(b, a, gl[a], gl[b], nt);
      wt[r] = nn[a] * nn[b] * dl * (dk - dl);
      ++r;
    }
  }
  kutil::ksum s;
  for (int q = 0; q < r; ++q) s.add(wt[q]);
  const double tot = s.val();
  for (int q = 0; q < r; ++q) wt[q] /= tot;
  return r;
}

int mtch(const double *f, int nu, int k, const int *trt, const int *ex, int ratio,
         double calip, int *mt, int *mc, double *md, double *wu, double *work) {
  // the covariance of the features, pooled within the two groups (as
  // MatchIt computes it)
  double *sg = work, *mn = work + z(k) * z(k), *dl = mn + 2 * z(k);
  int nn[2] = {0, 0};
  for (int i = 0; i < nu; ++i) ++nn[trt[i]];
  for (int grp = 0; grp < 2; ++grp) {
    for (int j = 0; j < k; ++j) {
      double s = 0.0;
      for (int i = 0; i < nu; ++i) s += trt[i] == grp ? f[z(j) * z(nu) + z(i)] : 0.0;
      mn[z(grp) * z(k) + z(j)] = s / nn[grp];
    }
  }
  for (std::size_t q = 0; q < z(k) * z(k); ++q) sg[q] = 0.0;
  for (int i = 0; i < nu; ++i) {
    const double *m = mn + z(trt[i]) * z(k);
    for (int j = 0; j < k; ++j) {
      for (int l = 0; l < k; ++l) {
        sg[z(l) * z(k) + z(j)] +=
            (f[z(j) * z(nu) + z(i)] - m[j]) * (f[z(l) * z(nu) + z(i)] - m[l]);
      }
    }
  }
  for (std::size_t q = 0; q < z(k) * z(k); ++q) sg[q] /= nu - 2.0;
  if (!chol(sg, k)) return -1;
  // each treated unit's ratio nearest comparison units in its exact group
  auto dist = [&](int i, int c) {
    // ||L^-1 (f_i - f_c)||
    double s2 = 0.0;
    for (int j = 0; j < k; ++j) {
      double s = f[z(j) * z(nu) + z(i)] - f[z(j) * z(nu) + z(c)];
      for (int l = 0; l < j; ++l) s -= sg[z(l) * z(k) + z(j)] * dl[l];
      dl[j] = s / sg[z(j) * z(k) + z(j)];
      s2 += dl[j] * dl[j];
    }
    return std::sqrt(s2);
  };
  int r = 0, nd = 0;
  kutil::ksum s1, s2;
  for (int i = 0; i < nu; ++i) {
    if (!trt[i]) continue;
    for (int c = 0; c < nu; ++c) {
      if (trt[c] || ex[c] != ex[i]) continue;
      const double d = dist(i, c);
      s1.add(d);
      s2.add(d * d);
      ++nd;
    }
    // the ratio nearest, one at a time (ties: the first unit)
    const int r0 = r;
    for (int q = 0; q < ratio; ++q) {
      int best = -1;
      double bd = 0.0;
      for (int c = 0; c < nu; ++c) {
        if (trt[c] || ex[c] != ex[i]) continue;
        bool taken = false;
        for (int a = r0; a < r; ++a) taken = taken || mc[a] == c;
        if (taken) continue;
        const double d = dist(i, c);
        if (best < 0 || d < bd) {
          best = c;
          bd = d;
        }
      }
      if (best < 0) break;
      mt[r] = i;
      mc[r] = best;
      md[r] = bd;
      ++r;
    }
  }
  if (r == 0) return 0;
  if (calip > 0.0) {
    const double mean = s1.val() / nd;
    const double sd = nd > 1 ? std::sqrt((s2.val() - nd * mean * mean) / (nd - 1)) : 0.0;
    int q = 0;
    for (int a = 0; a < r; ++a) {
      if (md[a] <= calip * sd) {
        mt[q] = mt[a];
        mc[q] = mc[a];
        md[q] = md[a];
        ++q;
      }
    }
    r = q;
  }
  for (int i = 0; i < nu; ++i) wu[i] = 0.0;
  for (int a = 0; a < r; ++a) {
    int km = 0;
    for (int b = 0; b < r; ++b) km += mt[b] == mt[a];
    wu[mt[a]] = 1.0;
    wu[mc[a]] += 1.0 / km;
  }
  return r;
}

}  // namespace pev

extern "C" int clab_didfe(const double *y, const double *x, int n, int p, const int *u, int nu,
                          const int *t, int nt, int fe, const double *w, const int *cl, int ncl,
                          const int *cl2, int ncl2, double tol, int maxit, double *beta,
                          double *vc, double *st, double *work, int *iwork) {
  if (n < 2 || p < 1 || nu < 1 || nt < 1 || ncl < 2 || (cl2 && ncl2 < 2)) return CLAB_EARG;
  const std::size_t sn = static_cast<std::size_t>(n), sp = static_cast<std::size_t>(p);
  double *v = work;              // y then the columns of x, demeaned: n (p + 1)
  double *r = v + sn * (sp + 1);  // residuals: n
  double *xi = r + sn;           // p p
  double *sc = xi + sp * sp;     // scratch
  for (std::size_t i = 0; i < sn; ++i) v[i] = y[i];
  for (std::size_t q = 0; q < sn * sp; ++q) v[sn + q] = x[q];
  const int sw = pev::within(v, n, p + 1, u, nu, t, nt, fe, w, tol, maxit, sc);
  if (sw < 0) return CLAB_EDOM;
  if (pev::ols(v + sn, v, n, p, w, beta, xi, r, sc) != CLAB_OK) return CLAB_EDOM;
  // fixest's "nested" count: an effect nested in a cluster costs nothing
  int *seen = iwork;
  auto nested = [&](const int *f, int nf, const int *c) {
    for (int a = 0; a < nf; ++a) seen[a] = -1;
    for (int i = 0; i < n; ++i) {
      if (seen[f[i]] < 0) seen[f[i]] = c[i];
      if (seen[f[i]] != c[i]) return false;
    }
    return true;
  };
  double kk = p + 1.0;
  if (fe && !(nested(u, nu, cl) || (cl2 && nested(u, nu, cl2)))) kk += nu - 1.0;
  if (!(nested(t, nt, cl) || (cl2 && nested(t, nt, cl2)))) kk += nt - 1.0;
  int *c12 = iwork, n12 = 0;
  if (cl2) n12 = pev::pairs(cl, cl2, n, ncl, ncl2, c12, iwork + sn);
  const int g = pev::crv(v + sn, r, w, n, p, xi, cl, ncl, cl2, ncl2, c12, n12, kk, vc, sc);
  kutil::ksum rss;
  for (int i = 0; i < n; ++i) rss.add((w ? w[i] : 1.0) * r[i] * r[i]);
  st[0] = sw;
  st[1] = g;
  st[2] = kk;
  st[3] = static_cast<double>(g) / (g - 1) * (n - 1.0) / (n - kk);
  st[4] = rss.val();
  return CLAB_OK;
}
