// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, forensic toxicology: calibration lines (clab_txcal),
// left-censored values (clab_lodimp), postmortem redistribution
// (clab_pmrrat), ethanol after death (clab_etohlr) and the reliability of
// specimen matrices (clab_mtxrel).
#include <algorithm>
#include <cmath>

#include "crimlab.h"
#include "kutil.h"

int clab_txcal(const double *x, const double *y, const double *w, int n, const double *u, int nu,
               double klod, double kloq, double *st, double *xh, int *flag) {
  if (n < 3 || nu < 0) return CLAB_EARG;
  if (!kutil::pos(klod) || !kutil::pos(kloq)) return CLAB_EDOM;
  kutil::ksum sw, sx, sy;
  for (int i = 0; i < n; ++i) {
    if (!kutil::pos(x[i]) || !std::isfinite(y[i]) || !kutil::pos(w[i])) return CLAB_EDOM;
    sw.add(w[i]);
    sx.add(w[i] * x[i]);
    sy.add(w[i] * y[i]);
  }
  const double xm = sx.val() / sw.val(), ym = sy.val() / sw.val();
  kutil::ksum sxx, sxy, syy;
  for (int i = 0; i < n; ++i) {
    sxx.add(w[i] * (x[i] - xm) * (x[i] - xm));
    sxy.add(w[i] * (x[i] - xm) * (y[i] - ym));
    syy.add(w[i] * (y[i] - ym) * (y[i] - ym));
  }
  if (!(sxx.val() > 0)) return CLAB_EDOM;
  const double b = sxy.val() / sxx.val(), a = ym - b * xm;
  if (b == 0) return CLAB_EDOM;
  kutil::ksum rss;
  for (int i = 0; i < n; ++i) {
    const double e = y[i] - a - b * x[i];
    rss.add(w[i] * e * e);
  }
  const double s = std::sqrt(rss.val() / (n - 2));
  st[0] = a;
  st[1] = b;
  st[2] = 1 - rss.val() / syy.val();  // syy > 0: the slope is not 0
  st[3] = s;
  st[4] = klod * s / std::fabs(b);
  st[5] = kloq * s / std::fabs(b);
  for (int j = 0; j < nu; ++j) {
    xh[j] = (u[j] - a) / b;
    flag[j] = std::isnan(xh[j]) ? -1 : xh[j] < st[4] ? 0 : xh[j] < st[5] ? 1 : 2;
  }
  return CLAB_OK;
}

int clab_lodimp(const double *x, int n, const double *lod, int nl, int meth, double *out,
                int *cens, double *st) {
  if (n < 1 || (nl != 1 && nl != n) || meth < 0 || meth > 2) return CLAB_EARG;
  for (int k = 0; k < nl; ++k) {
    if (!kutil::pos(lod[k])) return CLAB_EDOM;
  }
  const double f = meth == 0 ? 0.5 : meth == 1 ? 1 / std::sqrt(2.0) : 1.0;
  int meas = 0, nc = 0;
  for (int i = 0; i < n; ++i) {
    const double l = lod[nl == 1 ? 0 : i];
    cens[i] = !std::isnan(x[i]) && x[i] < l;
    out[i] = cens[i] ? f * l : x[i];
    meas += !std::isnan(x[i]);
    nc += cens[i];
  }
  st[0] = meas ? static_cast<double>(nc) / meas : kutil::nan();
  st[1] = n - meas;
  return CLAB_OK;
}

int clab_pmrrat(const double *c, const double *p, int n, double c1, double c2, double *r,
                int *cls) {
  if (n < 1) return CLAB_EARG;
  if (!kutil::pos(c1) || !std::isfinite(c2) || c2 < c1) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!kutil::pos(c[i]) || !kutil::pos(p[i])) return CLAB_EDOM;
    r[i] = c[i] / p[i];
    cls[i] = r[i] <= c1 ? 0 : r[i] <= c2 ? 1 : 2;
  }
  return CLAB_OK;
}

int clab_etohlr(const double *eth, const double *pr, const double *bu, const double *etg, int n,
                double thr, double *hi, int *v) {
  if (n < 1) return CLAB_EARG;
  if (!std::isfinite(thr)) return CLAB_EDOM;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(eth[i])) return CLAB_EDOM;
    const double a = pr[i], b = bu[i];
    hi[i] = std::isnan(a) ? b : std::isnan(b) ? a : std::max(a, b);
    v[i] = etg[i] > 0 ? 0 : hi[i] > thr ? 1 : 2;
  }
  return CLAB_OK;
}

int clab_mtxrel(const double *base, const int *dil, int n, double days, double stage,
                double dmax, double dcap, double smax, double scap, double *rel, int *ord) {
  if (n < 1) return CLAB_EARG;
  if (!kutil::nonneg(days) || !kutil::nonneg(stage) || !kutil::pos(dmax) || !kutil::pos(smax) ||
      !(dcap >= 0 && dcap <= 1) || !(scap >= 0 && scap <= 1)) {
    return CLAB_EDOM;
  }
  const double dp = 1 - std::min(days / dmax, dcap), sp = 1 - std::min(stage / smax, scap);
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(base[i])) return CLAB_EDOM;
    rel[i] = base[i] * sp * (dil[i] ? dp : 1.0);
    ord[i] = i;
  }
  // ponytail: insertion sort, stable; a body has ten specimen matrices
  for (int i = 1; i < n; ++i) {
    const int k = ord[i];
    int j = i;
    for (; j > 0 && rel[ord[j - 1]] < rel[k]; --j) ord[j] = ord[j - 1];
    ord[j] = k;
  }
  return CLAB_OK;
}
