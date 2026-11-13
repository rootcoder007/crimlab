// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, concentration of counts across units (clab_cconc), its
// Poisson null (clab_cnull), the growth of distinct places (clab_cgrow),
// the urban-scaling exponent (clab_urbsc), the Hurwitz zeta function
// (clab_hzta) and the tail index of a heavy-tailed count (clab_hilltl).
#include <cmath>
#include <cstdint>

#include "crimlab.h"
#include "kutil.h"

namespace {

// Gini coefficient of the units o[a..b) of x in ascending order, with
// weights w (NULL: 1 each): 1 - sum_i p_i (L_{i-1} + L_i), p the units'
// shares and L the cumulative shares of the total. 0 for no total.
double gini(const double *x, const double *w, const int *o, int a, int b) {
  kutil::ksum ws, ts;
  for (int i = a; i < b; ++i) {
    const double wi = w ? w[o[i]] : 1.0;
    ws.add(wi);
    ts.add(wi * x[o[i]]);
  }
  const double tw = ws.val(), tt = ts.val();
  if (!(tt > 0)) return 0.0;
  kutil::ksum s, run;
  double prev = 0.0;
  for (int i = a; i < b; ++i) {
    const double wi = w ? w[o[i]] : 1.0;
    run.add(wi * x[o[i]]);
    const double cur = run.val() / tt;
    s.add(wi / tw * (prev + cur));
    prev = cur;
  }
  return 1.0 - s.val();
}

// A Poisson draw: inversion below a mean of 10, otherwise the transformed
// rejection with squeeze of Hormann (1993, doi:10.1016/0167-6687(93)90997-4).
int rpois(kutil::rng &g, double mu) {
  if (mu < 10) {
    const double e = std::exp(-mu);
    double p = e, c = e, u = g.unif();
    int k = 0;
    while (u > c && k < 1000) {
      ++k;
      p *= mu / k;
      c += p;
    }
    return k;
  }
  const double sl = std::sqrt(mu), ll = std::log(mu);
  const double b = 0.931 + 2.53 * sl, a = -0.059 + 0.02483 * b;
  const double ia = 1.1239 + 1.1328 / (b - 3.4), vr = 0.9277 - 3.6224 / (b - 2);
  for (;;) {
    const double u = g.unif() - 0.5, v = g.unif(), us = 0.5 - std::fabs(u);
    const double k = std::floor((2 * a / us + b) * u + mu + 0.43);
    if (us >= 0.07 && v <= vr) return static_cast<int>(k);
    if (k < 0 || (us < 0.013 && v > us)) continue;
    if (std::log(v) + std::log(ia) - std::log(a / (us * us) + b) <=
        -mu + k * ll - std::lgamma(k + 1)) {
      return static_cast<int>(k);
    }
  }
}

// The Bernoulli numbers B_2, B_4, ..., B_20 over (2j)!.
const double bern[10] = {1.0 / 6 / 2,
                         -1.0 / 30 / 24,
                         1.0 / 42 / 720,
                         -1.0 / 30 / 40320,
                         5.0 / 66 / 3628800,
                         -691.0 / 2730 / 479001600,
                         7.0 / 6 / 87178291200.0,
                         -3617.0 / 510 / 20922789888000.0,
                         43867.0 / 798 / 6402373705728000.0,
                         -174611.0 / 330 / 2432902008176640000.0};

// Hurwitz zeta, s > 1, q > 0: the first N terms summed, then the
// Euler-Maclaurin tail at a = q + N >= max(20, 2 s), ten Bernoulli terms.
double hz(double s, double q) {
  const double lim = std::fmax(20.0, 2.0 * s);
  const int nn = q >= lim ? 0 : static_cast<int>(std::ceil(lim - q));
  kutil::ksum t;
  for (int k = nn - 1; k >= 0; --k) t.add(std::pow(q + k, -s));
  const double a = q + nn;
  const double ap = std::pow(a, -s);
  t.add(a * ap / (s - 1));
  t.add(ap / 2);
  // s (s + 1) ... (s + 2j - 2) a^(-s - 2j + 1)
  double f = s * ap / a;
  for (int j = 0; j < 10; ++j) {
    t.add(bern[j] * f);
    f *= (s + 2 * j + 1) * (s + 2 * j + 2) / (a * a);
  }
  return t.val();
}

// The negative log-likelihood of the discrete power law on the tail, up
// to a constant: n log zeta(a, xmin) + a sum log x.
double nll(double a, double n, double ls, double xmin) {
  return n * std::log(hz(a, xmin)) + a * ls;
}

}  // namespace

int clab_cconc(const double *x, const double *w, int n, const double *fr, int nf, double *st,
               double *lp, double *lv, double *top, double *topn, int *o) {
  if (n < 1 || nf < 0) return CLAB_EARG;
  kutil::ksum ws, ts, zs;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(x[i]) || (w && !kutil::pos(w[i]))) return CLAB_EDOM;
    const double wi = w ? w[i] : 1.0;
    ws.add(wi);
    ts.add(wi * x[i]);
    if (x[i] == 0) zs.add(wi);
    o[i] = i;
  }
  for (int f = 0; f < nf; ++f) {
    if (!(fr[f] >= 0 && fr[f] <= 1)) return CLAB_EDOM;
  }
  const double tw = ws.val(), tt = ts.val();
  if (!(tt > 0)) return CLAB_EDOM;
  kutil::hsort(o, n, [&](int a, int b) { return x[a] < x[b] || (x[a] == x[b] && a < b); });
  const double m = tt / tw;
  kutil::ksum ss;
  for (int i = 0; i < n; ++i) ss.add((w ? w[i] : 1.0) * (x[i] - m) * (x[i] - m));
  const double v = tw > 1 ? ss.val() / (tw - 1) : kutil::nan(), vp = ss.val() / tw;
  int z0 = 0;
  while (z0 < n && x[o[z0]] == 0) ++z0;
  const double zs0 = zs.val() / tw, zn = std::exp(-m);
  const double g = gini(x, w, o, 0, n), gp = gini(x, w, o, z0, n);
  st[0] = tw;
  st[1] = tt;
  st[2] = m;
  st[3] = v;
  st[4] = v / m;
  st[5] = std::fmax(v - m, 0.0);
  st[6] = std::sqrt(st[5]);
  st[7] = zs0;
  st[8] = zn;
  st[9] = zs0 - zn;
  st[10] = g;
  st[11] = gp;
  st[12] = zn + (1 - zn) * gp;
  st[13] = vp;
  st[14] = std::sqrt(vp / 2) / m;
  st[15] = 1 / std::sqrt(2 * m);
  // the Lorenz curve, one point per unit in ascending order
  lp[0] = 0.0;
  lv[0] = 0.0;
  kutil::ksum pc, vc;
  for (int i = 0; i < n; ++i) {
    const double wi = w ? w[o[i]] : 1.0;
    pc.add(wi);
    vc.add(wi * x[o[i]]);
    lp[i + 1] = pc.val() / tw;
    lv[i + 1] = vc.val() / tt;
  }
  // top shares: the largest ceil(f W) units, a unit of weight w counting
  // as w units (the last one in part)
  for (int f = 0; f < nf; ++f) {
    const double k = std::fmin(tw, std::ceil(fr[f] * tw - 1e-9));
    double left = k;
    kutil::ksum s;
    for (int i = n - 1; i >= 0 && left > 0; --i) {
      const double take = std::fmin(left, w ? w[o[i]] : 1.0);
      s.add(take * x[o[i]]);
      left -= take;
    }
    top[f] = s.val() / tt;
    topn[f] = k;
  }
  return CLAB_OK;
}

int clab_cnull(double mu, int n, int nsim, double seed, double *g, int *work) {
  if (n < 1 || nsim < 0) return CLAB_EARG;
  if (!kutil::pos(mu) || mu > 1e8 || !(seed >= 0)) return CLAB_EDOM;
  kutil::rng r(static_cast<std::uint64_t>(seed));
  for (int s = 0; s < nsim; ++s) {
    for (int i = 0; i < n; ++i) work[i] = rpois(r, mu);
    kutil::hsort(work, n, [](int a, int b) { return a < b; });
    // 2 sum i x_(i) / (n T) - (n + 1) / n, 0 for no events
    double t = 0, c = 0;
    for (int i = 0; i < n; ++i) {
      t += work[i];
      c += static_cast<double>(i + 1) * work[i];
    }
    g[s] = t > 0 ? 2 * c / (static_cast<double>(n) * t) - (n + 1.0) / n : 0.0;
  }
  return CLAB_OK;
}

int clab_cgrow(const int *code, int n, int ng, double *dist, double *mhat, double *ex, double *lo,
               double *hi, double *slope, int *seen) {
  if (n < 2 || ng < 1) return CLAB_EARG;
  for (int j = 0; j < ng; ++j) seen[j] = 0;
  int k = 0;
  for (int i = 0; i < n; ++i) {
    if (code[i] < 0 || code[i] >= ng) return CLAB_EARG;
    if (!seen[code[i]]++) ++k;
    dist[i] = k;
  }
  // M solves M sum_{j < n} 1 / (M + j) = K, which rises from 1 to n in M
  auto f = [&](double mm) {
    kutil::ksum s;
    for (int j = n - 1; j >= 0; --j) s.add(mm / (mm + j));
    return s.val();
  };
  double m;
  if (k >= n) {
    m = INFINITY;
  } else if (k <= 1) {
    m = 0.0;
  } else {
    double a = 1e-8, b = 1e8;
    while (f(b) < k) b *= 16;
    for (int it = 0; it < 200 && b / a - 1 > 1e-15; ++it) {
      const double c = std::sqrt(a * b);
      if (f(c) < k)
        a = c;
      else
        b = c;
    }
    m = std::sqrt(a * b);
  }
  *mhat = m;
  const bool fin = std::isfinite(m) && m > 0;
  kutil::ksum s;
  for (int i = 0; i < n; ++i) {
    s.add(fin ? 1 / (m + i) : 0.0);
    ex[i] = fin ? m * s.val() : kutil::nan();
    lo[i] = fin ? m * std::log1p((i + 1) / m) : kutil::nan();
    hi[i] = fin ? 1 + m * std::log1p(i / m) : kutil::nan();
  }
  // the slope of log distinct on log index over the second half
  kutil::ksum sx, sy;
  int c = 0;
  for (int i = 0; i < n; ++i) {
    if (2 * (i + 1) < n) continue;
    sx.add(std::log(i + 1.0));
    sy.add(std::log(dist[i]));
    ++c;
  }
  if (c > 2) {
    const double mx = sx.val() / c, my = sy.val() / c;
    kutil::ksum sxy, sxx;
    for (int i = 0; i < n; ++i) {
      if (2 * (i + 1) < n) continue;
      const double dx = std::log(i + 1.0) - mx;
      sxy.add(dx * (std::log(dist[i]) - my));
      sxx.add(dx * dx);
    }
    *slope = sxy.val() / sxx.val();
  } else {
    *slope = kutil::nan();
  }
  return CLAB_OK;
}

int clab_urbsc(const double *x, const double *pop, int n, double *out) {
  if (n < 1) return CLAB_EARG;
  kutil::ksum sx, sy;
  int c = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::nonneg(x[i]) || !kutil::pos(pop[i])) return CLAB_EDOM;
    if (x[i] == 0) continue;
    sx.add(std::log(pop[i]));
    sy.add(std::log(x[i]));
    ++c;
  }
  out[4] = c;
  if (c < 3) return CLAB_EDOM;
  const double mx = sx.val() / c, my = sy.val() / c;
  kutil::ksum sxy, sxx, syy;
  for (int i = 0; i < n; ++i) {
    if (x[i] == 0) continue;
    const double dx = std::log(pop[i]) - mx, dy = std::log(x[i]) - my;
    sxy.add(dx * dy);
    sxx.add(dx * dx);
    syy.add(dy * dy);
  }
  if (!(sxx.val() > 0)) return CLAB_EDOM;
  const double b = sxy.val() / sxx.val();
  kutil::ksum rs;
  for (int i = 0; i < n; ++i) {
    if (x[i] == 0) continue;
    const double r = std::log(x[i]) - my - b * (std::log(pop[i]) - mx);
    rs.add(r * r);
  }
  out[0] = b;
  out[1] = std::exp(my - b * mx);
  out[2] = std::sqrt(rs.val() / (c - 2) / sxx.val());
  out[3] = syy.val() > 0 ? 1 - rs.val() / syy.val() : kutil::nan();
  return CLAB_OK;
}

int clab_hzta(double s, double q, double *z) {
  if (!(s > 1) || !std::isfinite(s) || !kutil::pos(q)) return CLAB_EDOM;
  *z = hz(s, q);
  return CLAB_OK;
}

int clab_hilltl(const double *x, int n, double xmin, int discrete, int approx, double *out,
                double *work) {
  if (n < 1) return CLAB_EARG;
  if (!kutil::pos(xmin)) return CLAB_EDOM;
  int m = 0;
  for (int i = 0; i < n; ++i) {
    if (!kutil::pos(x[i])) return CLAB_EDOM;
    if (x[i] >= xmin) work[m++] = x[i];
  }
  out[3] = m;
  out[4] = 0.0;
  if (m < 2) return CLAB_EDOM;
  kutil::hsort(work, m, [](double a, double b) { return a < b; });
  if (work[0] == work[m - 1]) return CLAB_EDOM;
  double alpha, se;
  // the fitted distribution function at u, F(u) = P(X <= u)
  auto ks = [&](auto cdf, bool cont) {
    double d = 0;
    for (int i = 0; i < m;) {
      int j = i;
      while (j < m && work[j] == work[i]) ++j;
      const double fu = cdf(work[i]);
      d = std::fmax(d, std::fabs(static_cast<double>(j) / m - fu));
      if (cont) d = std::fmax(d, std::fabs(static_cast<double>(i) / m - fu));
      i = j;
    }
    return d;
  };
  if (discrete && !approx) {
    kutil::ksum l;
    for (int i = 0; i < m; ++i) l.add(std::log(work[i]));
    const double ls = l.val();
    // nll is convex in alpha (log zeta is a log-sum-exp of lines in
    // alpha): golden-section search on [1.0001, 25]
    const double r = (std::sqrt(5.0) - 1) / 2;
    double a = 1.0001, b = 25, c = b - r * (b - a), d = a + r * (b - a);
    double fc = nll(c, m, ls, xmin), fd = nll(d, m, ls, xmin);
    while (b - a > 1e-11 * (1 + std::fabs(c))) {
      if (fc < fd) {
        b = d;
        d = c;
        fd = fc;
        c = b - r * (b - a);
        fc = nll(c, m, ls, xmin);
      } else {
        a = c;
        c = d;
        fc = fd;
        d = a + r * (b - a);
        fd = nll(d, m, ls, xmin);
      }
    }
    alpha = (a + b) / 2;
    out[4] = alpha <= 1.0001 + 1e-3 || alpha >= 25 - 1e-3 ? 1.0 : 0.0;
    // the observed information by a central difference, step 1e-4
    const double h = 1e-4;
    const double d2 =
        (nll(alpha + h, m, ls, xmin) - 2 * nll(alpha, m, ls, xmin) + nll(alpha - h, m, ls, xmin)) /
        (h * h);
    se = d2 > 0 ? std::sqrt(1 / d2) : kutil::nan();
    const double zq = hz(alpha, xmin);
    out[2] = ks([&](double u) { return 1 - hz(alpha, u + 1) / zq; }, false);
  } else {
    // the continuous Hill estimator; for counts, with xmin - 1/2
    const double den = discrete && xmin > 0.5 ? xmin - 0.5 : xmin;
    kutil::ksum l;
    for (int i = 0; i < m; ++i) l.add(std::log(work[i] / den));
    alpha = 1 + m / l.val();
    se = (alpha - 1) / std::sqrt(static_cast<double>(m));
    out[2] = ks([&](double u) { return 1 - std::pow(u / den, 1 - alpha); }, !discrete);
  }
  out[0] = alpha;
  out[1] = se;
  return CLAB_OK;
}
