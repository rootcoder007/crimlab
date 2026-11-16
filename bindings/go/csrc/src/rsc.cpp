// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// The Rao-Scott test of independence in a two-way table from a complex
// survey (crimlab.h, rscchi()): Pearson's X^2 on the weighted cell
// proportions scaled to the sample size, and its corrections by the
// generalised design effects delta, the eigenvalues of Delta = (C' D^-1
// C / n)^-1 (C' D^-1 V D^-1 C), C a basis of the interaction contrasts,
// D the diagonal of the cell proportions and V their design covariance
// by linearisation over the first-stage units within strata (Rao and
// Scott 1981, 1984; the estimator of R's survey package, svychisq()).
// First order: X^2 / mean(delta) on (r - 1)(c - 1) degrees of freedom;
// second order: X^2 / (mean(delta) (1 + a^2)), a the coefficient of
// variation of delta, on (r - 1)(c - 1) / (1 + a^2), the Satterthwaite
// form; and the F form, X^2 / tr(Delta) on d0 and d0 nu, d0 = tr(Delta)^2
// / tr(Delta^2), nu the first-stage units less the strata.
#include <cmath>
#include <cstddef>
#include <vector>

#include "crimlab.h"
#include "gla.h"
#include "kutil.h"

extern "C" int clab_rscchi(const int *r, const int *c, int nr, int nc, int n, const double *w, const int *st,
                           int nst, const int *ps, int nps, const double *fpc, double *tab, double *test,
                           double *dlt, double *aux) CLAB_TRY {
  using std::size_t;
  auto z = [](int i) { return static_cast<size_t>(i); };
  if (nr < 2 || nc < 2 || n < 2 || nst < 1 || nps <= nst) return CLAB_EARG;
  const int K0 = nr * nc, K = (nr - 1) * (nc - 1);
  // the cell proportions, and the weights of each first-stage unit by
  // cell and in all
  std::vector<double> a(z(nps) * z(K0), 0.0), wa(z(nps), 0.0), pk(z(K0), 0.0);
  std::vector<int> sp(z(nps), -1);
  double sw = 0;
  for (int i = 0; i < n; ++i) {
    if (r[i] < 0 || r[i] >= nr || c[i] < 0 || c[i] >= nc || ps[i] < 0 || ps[i] >= nps || st[i] < 0 ||
        st[i] >= nst || !(w[i] > 0)) {
      return CLAB_EARG;
    }
    const int k = r[i] + nr * c[i];
    pk[z(k)] += w[i];
    a[z(ps[i]) * z(K0) + z(k)] += w[i];
    wa[z(ps[i])] += w[i];
    sp[z(ps[i])] = st[i];
    sw += w[i];
  }
  for (double &v : pk) v /= sw;
  std::vector<double> pr(z(nr), 0.0), pc(z(nc), 0.0);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      pr[z(i)] += pk[z(i + nr * j)];
      pc[z(j)] += pk[z(i + nr * j)];
    }
  }
  double x2 = 0;
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      const double e = pr[z(i)] * pc[z(j)], dv = pk[z(i + nr * j)] - e;
      x2 += dv * dv / e;
    }
  }
  x2 *= n;
  for (int k = 0; k < K0; ++k) tab[k] = pk[z(k)];
  // V: the covariance of the cell proportions, from the first-stage
  // totals of the influence values w (I - p) / sum w
  std::vector<double> nh(z(nst), 0.0);
  for (int q = 0; q < nps; ++q) nh[z(sp[z(q)])] += 1;
  for (double v : nh) {
    if (v < 2) return CLAB_EDOM;
  }
  std::vector<double> t(z(nps) * z(K0)), mh(z(nst) * z(K0), 0.0);
  for (int q = 0; q < nps; ++q) {
    for (int k = 0; k < K0; ++k) {
      const double v = (a[z(q) * z(K0) + z(k)] - pk[z(k)] * wa[z(q)]) / sw;
      t[z(q) * z(K0) + z(k)] = v;
      mh[z(sp[z(q)]) * z(K0) + z(k)] += v / nh[z(sp[z(q)])];
    }
  }
  std::vector<double> V(z(K0) * z(K0), 0.0);
  for (int q = 0; q < nps; ++q) {
    if (kutil::stop(static_cast<double>(K0) * K0)) return CLAB_EINTR;
    const int h = sp[z(q)];
    const double f = fpc ? nh[z(h)] / fpc[h] : 0.0, m = (1 - f) * nh[z(h)] / (nh[z(h)] - 1);
    for (int k = 0; k < K0; ++k) {
      const double dk = t[z(q) * z(K0) + z(k)] - mh[z(h) * z(K0) + z(k)];
      for (int l = 0; l <= k; ++l) {
        V[z(k) + z(l) * z(K0)] += m * dk * (t[z(q) * z(K0) + z(l)] - mh[z(h) * z(K0) + z(l)]);
      }
    }
  }
  // the contrasts (e_i - 1 / nr) x (f_j - 1 / nc), i < nr - 1, j < nc - 1:
  // a basis of the interaction space, the space the main effects leave
  std::vector<double> C(z(K0) * z(K)), id(z(K0));
  for (int i = 0; i + 1 < nr; ++i) {
    for (int j = 0; j + 1 < nc; ++j) {
      const int col = i + (nr - 1) * j;
      for (int rr = 0; rr < nr; ++rr) {
        for (int cc = 0; cc < nc; ++cc) {
          C[z(col) * z(K0) + z(rr + nr * cc)] = ((rr == i ? 1.0 : 0.0) - 1.0 / nr) * ((cc == j ? 1.0 : 0.0) - 1.0 / nc);
        }
      }
    }
  }
  for (int k = 0; k < K0; ++k) id[z(k)] = pk[z(k)] > 0 ? 1 / pk[z(k)] : 0.0;
  // den = C' D^-1 C / n; num = C' D^-1 V D^-1 C
  std::vector<double> dc(z(K0) * z(K)), vdc(z(K0) * z(K)), den(z(K) * z(K)), num(z(K) * z(K));
  for (int u = 0; u < K; ++u) {
    for (int k = 0; k < K0; ++k) dc[z(u) * z(K0) + z(k)] = id[z(k)] * C[z(u) * z(K0) + z(k)];
  }
  for (int u = 0; u < K; ++u) {
    for (int k = 0; k < K0; ++k) {
      double s = 0;
      for (int l = 0; l < K0; ++l) {
        const double v = l <= k ? V[z(k) + z(l) * z(K0)] : V[z(l) + z(k) * z(K0)];
        s += v * dc[z(u) * z(K0) + z(l)];
      }
      vdc[z(u) * z(K0) + z(k)] = s;
    }
  }
  for (int u = 0; u < K; ++u) {
    for (int v = 0; v < K; ++v) {
      double s1 = 0, s2 = 0;
      for (int k = 0; k < K0; ++k) {
        s1 += C[z(u) * z(K0) + z(k)] * dc[z(v) * z(K0) + z(k)];
        s2 += dc[z(u) * z(K0) + z(k)] * vdc[z(v) * z(K0) + z(k)];
      }
      den[z(u) + z(v) * z(K)] = s1 / n;
      num[z(u) + z(v) * z(K)] = s2;
    }
  }
  // Delta's eigenvalues: those of L^-1 num L^-T, den = L L'
  if (!gla::chol(den.data(), K)) return CLAB_EDOM;
  gla::fsolve(den.data(), K, num.data(), K);
  std::vector<double> m(z(K) * z(K));
  for (int u = 0; u < K; ++u) {
    for (int v = 0; v < K; ++v) m[z(v) + z(u) * z(K)] = num[z(u) + z(v) * z(K)];
  }
  gla::fsolve(den.data(), K, m.data(), K);
  for (int u = 0; u < K; ++u) {
    for (int v = 0; v < u; ++v) {
      const double s = (m[z(u) + z(v) * z(K)] + m[z(v) + z(u) * z(K)]) / 2;
      m[z(u) + z(v) * z(K)] = m[z(v) + z(u) * z(K)] = s;
    }
  }
  std::vector<double> vec(z(K) * z(K));
  if (!gla::jacobi(m.data(), K, dlt, vec.data())) return CLAB_EDOM;
  double tr = 0, tr2 = 0;
  for (int u = 0; u < K; ++u) {
    tr += dlt[u];
    tr2 += dlt[u] * dlt[u];
  }
  const double db = tr / K, a2 = (tr2 / K - db * db) / (db * db), d0 = tr * tr / tr2, nu = nps - nst;
  if (!(db > 0)) return CLAB_EDOM;
  // pearson, first order, second order, F: stat, df1, df2, p
  const double s[4] = {x2, x2 / db, x2 / (db * (1 + a2)), x2 / tr};
  const double f1[4] = {static_cast<double>(K), static_cast<double>(K), K / (1 + a2), d0};
  for (int k = 0; k < 4; ++k) {
    test[k] = s[k];
    test[4 + k] = f1[k];
    test[8 + k] = k == 3 ? d0 * nu : kutil::nan();
    test[12 + k] = k == 3 ? clab_pf(s[k], d0, d0 * nu, 0) : clab_pchisq(s[k], f1[k], 0);
  }
  aux[0] = db;
  aux[1] = a2;
  aux[2] = nu;
  return CLAB_OK;
}
CLAB_CATCH(CLAB_ENOMEM)
