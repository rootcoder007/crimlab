// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, contingency tables: the Pearson and likelihood-ratio
// statistics with the association measures (clab_ctabx), and Fisher's
// exact test with, for a 2 x 2 table, the conditional maximum likelihood
// odds ratio and its exact interval (clab_fshx). No allocation.
#include <cfloat>
#include <cmath>
#include <limits>

#include "crimlab.h"

namespace {

constexpr double qnan = std::numeric_limits<double>::quiet_NaN();
constexpr double inf = std::numeric_limits<double>::infinity();

double at(const double *x, int r, int i, int j) {
  return x[static_cast<long long>(j) * r + i];
}

// ------------------------------------------------- Fisher, r x c

// The tables with the margins of the observed one, cell by cell in
// column-major order; each adds its probability when it is at most the
// observed one's (to a relative 1e-7, as R).
struct fen {
  int r, c;
  double *rem;     // r: what is left of each row margin
  const double *col;
  double lcut, lconst, p, cnt, cap;
  bool over;
};

// cell (i, j) with s left of column j's margin; lp: the sum of
// lgamma(cell + 1) so far
void walk(fen &f, int i, int j, double s, double lp) {
  if (f.over) return;
  if (j == f.c - 1) {
    // the last column is what the rows have left
    double l = lp;
    for (int q = 0; q < f.r; ++q) l += std::lgamma(f.rem[q] + 1);
    if (++f.cnt > f.cap) {
      f.over = true;
      return;
    }
    const double lt = f.lconst - l;
    if (lt <= f.lcut) f.p += std::exp(lt);
    return;
  }
  if (i == f.r - 1) {
    // the last cell of the column is what the column has left
    if (s > f.rem[i]) return;
    f.rem[i] -= s;
    walk(f, 0, j + 1, f.col[j + 1], lp + std::lgamma(s + 1));
    f.rem[i] += s;
    return;
  }
  // the rows below must be able to take the rest
  double below = 0;
  for (int q = i + 1; q < f.r; ++q) below += f.rem[q];
  const double lo = std::fmax(0.0, s - below), hi = std::fmin(s, f.rem[i]);
  for (double v = lo; v <= hi; v += 1) {
    f.rem[i] -= v;
    walk(f, i + 1, j, s - v, lp + std::lgamma(v + 1));
    f.rem[i] += v;
    if (f.over) return;
  }
}

// ------------------------------------------------- Fisher, 2 x 2

// The noncentral hypergeometric distribution of the (1, 1) cell given
// the margins (m in the first column, nn in the second, k in the first
// row) and the odds ratio t: R's dnhyper, mnhyper and pnhyper.
struct nh {
  double m, nn, k, lo, hi;
  double ld(double x) const {
    return std::lgamma(m + 1) - std::lgamma(x + 1) - std::lgamma(m - x + 1) +
           std::lgamma(nn + 1) - std::lgamma(k - x + 1) - std::lgamma(nn - k + x + 1);
  }
  // sum over the support of exp(ld(x) + x log t - max) f(x), normalised
  template <class F>
  double ex(double t, F f) const {
    double mx = -inf;
    for (double x = lo; x <= hi; x += 1) mx = std::fmax(mx, ld(x) + x * std::log(t));
    double s = 0, a = 0;
    for (double x = lo; x <= hi; x += 1) {
      const double d = std::exp(ld(x) + x * std::log(t) - mx);
      s += d;
      a += d * f(x);
    }
    return a / s;
  }
  double mean(double t) const {
    if (t == 0) return lo;
    if (std::isinf(t)) return hi;
    return ex(t, [](double x) { return x; });
  }
  double cdf(double q, double t, bool up) const {
    if (t == 0) return up ? (q <= lo ? 1.0 : 0.0) : 1.0;
    if (std::isinf(t)) return up ? 1.0 : (q >= hi ? 1.0 : 0.0);
    return ex(t, [q, up](double x) { return (up ? x >= q : x <= q) ? 1.0 : 0.0; });
  }
};

// a root of the increasing (inc) or decreasing g on (0, 1) by bisection
template <class G>
double root01(G g, bool inc) {
  double a = 0, b = 1;
  for (int it = 0; it < 200 && b - a > 4 * DBL_EPSILON * b; ++it) {
    const double mid = (a + b) / 2;
    if ((g(mid) < 0) == inc) a = mid;
    else b = mid;
  }
  return (a + b) / 2;
}

}  // namespace

extern "C" {

int clab_ctabx(const double *x, int r, int c, int yates, const int *dg, double *e,
               double *res, double *st) {
  if (r < 1 || c < 1) return CLAB_EARG;
  const long long rc = static_cast<long long>(r) * c;
  double n = 0;
  for (long long q = 0; q < rc; ++q) {
    if (!(x[q] >= 0) || !std::isfinite(x[q])) return CLAB_EDOM;
    n += x[q];
  }
  if (!(n > 0)) return CLAB_EDOM;
  double ch = 0, chy = 0, g2 = 0, d = 0;
  bool anyd = false;
  for (int i = 0; i < r; ++i) {
    double ri = 0;
    for (int j = 0; j < c; ++j) ri += at(x, r, i, j);
    if (!(ri > 0)) return CLAB_EDOM;
    for (int j = 0; j < c; ++j) {
      double cj = 0;
      for (int q = 0; q < r; ++q) cj += at(x, r, q, j);
      if (!(cj > 0)) return CLAB_EDOM;
      const long long k = static_cast<long long>(j) * r + i;
      const double o = x[k], ex = ri * cj / n;
      e[k] = ex;
      res[k] = (o - ex) / std::sqrt(ex);
      ch += (o - ex) * (o - ex) / ex;
      const double y = std::fabs(o - ex) - std::fmin(0.5, std::fabs(o - ex));
      chy += y * y / ex;
      if (o > 0) g2 += 2 * o * std::log(o / ex);
    }
    if (dg && dg[i] >= 0) {
      if (dg[i] >= c) return CLAB_EARG;
      d += at(x, r, i, dg[i]);
      anyd = true;
    }
  }
  const int k = r < c ? r : c;
  const double df = static_cast<double>(r - 1) * (c - 1);
  st[0] = n;
  st[1] = yates && r == 2 && c == 2 ? chy : ch;
  st[2] = df;
  st[3] = g2;
  st[4] = k > 1 ? std::sqrt(ch / (n * (k - 1))) : qnan;
  // Bergsma's (2013) bias correction of phi^2 and of the dimensions
  if (k > 1 && n > 1) {
    const double p2 = std::fmax(0.0, ch / n - df / (n - 1));
    const double rr = r - (r - 1.0) * (r - 1.0) / (n - 1), cc = c - (c - 1.0) * (c - 1.0) / (n - 1);
    const double den = std::fmin(rr - 1, cc - 1);
    st[5] = den > 0 ? std::sqrt(p2 / den) : qnan;
  } else {
    st[5] = qnan;
  }
  st[6] = std::sqrt(ch / (ch + n));
  if (r == 2 && c == 2) {
    const double a = x[0], b = x[2], cq = x[1], dd = x[3];
    st[7] = (a * dd - b * cq) /
            std::sqrt((a + b) * (cq + dd) * (a + cq) * (b + dd));
  } else {
    st[7] = qnan;
  }
  st[8] = anyd ? d / n : qnan;
  return CLAB_OK;
}

int clab_fshx(const double *x, int r, int c, double cap, double conf, double *st,
              double *rem) {
  if (r < 2 || c < 2 || !(cap >= 1) || !(conf > 0 && conf < 1)) return CLAB_EARG;
  const long long rc = static_cast<long long>(r) * c;
  double n = 0, lobs = 0;
  for (long long q = 0; q < rc; ++q) {
    if (!(x[q] >= 0) || x[q] != std::floor(x[q]) || !std::isfinite(x[q])) return CLAB_EDOM;
    n += x[q];
    lobs += std::lgamma(x[q] + 1);
  }
  double lc = -std::lgamma(n + 1);
  for (int i = 0; i < r; ++i) {
    double s = 0;
    for (int j = 0; j < c; ++j) s += at(x, r, i, j);
    rem[i] = s;
    lc += std::lgamma(s + 1);
  }
  double *col = rem + r;
  for (int j = 0; j < c; ++j) {
    double s = 0;
    for (int i = 0; i < r; ++i) s += at(x, r, i, j);
    col[j] = s;
    lc += std::lgamma(s + 1);
  }
  fen f{r, c, rem, col, lc - lobs + std::log1p(1e-7), lc, 0, 0, cap, false};
  walk(f, 0, 0, col[0], 0);
  if (f.over) return CLAB_EDOM;
  st[0] = std::fmin(1.0, f.p);
  st[4] = f.cnt;
  st[1] = st[2] = st[3] = qnan;
  if (r != 2 || c != 2) return CLAB_OK;
  // the 2 x 2 table: fisher.test()'s estimate and interval
  const double a = x[0];
  nh h{x[0] + x[1], x[2] + x[3], x[0] + x[2], 0, 0};
  h.lo = std::fmax(0.0, h.k - h.nn);
  h.hi = std::fmin(h.k, h.m);
  const double al = (1 - conf) / 2;
  // the conditional maximum likelihood estimate: the mean equals a
  if (a == h.lo) {
    st[1] = 0;
  } else if (a == h.hi) {
    st[1] = inf;
  } else {
    const double mu = h.mean(1);
    if (mu > a) st[1] = root01([&](double t) { return h.mean(t) - a; }, true);
    else st[1] = mu < a ? 1 / root01([&](double t) { return h.mean(1 / t) - a; }, false) : 1;
  }
  // lower end
  if (a == h.lo) {
    st[2] = 0;
  } else {
    const double p = h.cdf(a, 1, true);
    if (p > al) st[2] = root01([&](double t) { return h.cdf(a, t, true) - al; }, true);
    else st[2] = p < al ? 1 / root01([&](double t) { return h.cdf(a, 1 / t, true) - al; }, false) : 1;
  }
  // upper end
  if (a == h.hi) {
    st[3] = inf;
  } else {
    const double p = h.cdf(a, 1, false);
    if (p < al) st[3] = root01([&](double t) { return h.cdf(a, t, false) - al; }, false);
    else st[3] = p > al ? 1 / root01([&](double t) { return h.cdf(a, 1 / t, false) - al; }, true) : 1;
  }
  return CLAB_OK;
}

}  // extern "C"
