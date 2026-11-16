// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// Data description: histogram bins as R's hist() counts them, pairwise
// rank and product-moment correlations, runs and tallies of missing
// values, column moments, quantiles and robust spread, Mahalanobis
// distances with a classical or a median and MAD centre and scale,
// Little's test of missing completely at random by EM, and Cochran's Q.
// The sums run as the reference implementation and R run them: long
// double where R's mean(), colMeans() and cov() accumulate in long double,
// double where its kernels do.
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "crimlab.h"
#include "kutil.h"

namespace {

typedef long double ld;
using std::size_t;
using std::vector;

// R's mean(): the long double sum, refined by a second pass.
double rmean(const double *x, size_t n) {
  if (n == 0) return NAN;
  ld s = 0;
  for (size_t i = 0; i < n; ++i) s += x[i];
  s /= n;
  if (std::isfinite(static_cast<double>(s))) {
    ld t = 0;
    for (size_t i = 0; i < n; ++i) t += x[i] - s;
    s += t / n;
  }
  return static_cast<double>(s);
}

// The type-7 quantile of the reference kernel on sorted v: v[lo] + g (v[lo
// + 1] - v[lo]), the weighted mean when that difference is not finite.
double q7s(const vector<double> &v, double p) {
  const size_t n = v.size();
  if (n == 0) return NAN;
  if (std::isnan(p)) return p;
  if (n == 1) return v[0];
  const double h = (static_cast<double>(n) - 1.0) * p, fl = std::floor(h);
  long lo = static_cast<long>(fl);
  if (lo < 0) lo = 0;
  if (lo >= static_cast<long>(n) - 1) return v[n - 1];
  const double g = h - fl;
  if (g == 0.0) return v[static_cast<size_t>(lo)];
  const double d = v[static_cast<size_t>(lo) + 1] - v[static_cast<size_t>(lo)];
  if (std::isfinite(d)) return v[static_cast<size_t>(lo)] + g * d;
  return (1.0 - g) * v[static_cast<size_t>(lo)] + g * v[static_cast<size_t>(lo) + 1];
}

double med(const double *x, size_t n) {
  vector<double> v(x, x + n);
  std::sort(v.begin(), v.end());
  return q7s(v, 0.5);
}

double mad(const double *x, size_t n) {
  vector<double> v(x, x + n);
  std::sort(v.begin(), v.end());
  const double m = q7s(v, 0.5);
  for (size_t i = 0; i < n; ++i) v[i] = std::fabs(x[i] - m);
  std::sort(v.begin(), v.end());
  return 1.4826 * q7s(v, 0.5);
}

// Midranks, NaN last.
void midrank(const double *a, size_t n, double *out) {
  vector<size_t> o(n);
  for (size_t i = 0; i < n; ++i) o[i] = i;
  std::sort(o.begin(), o.end(), [&](size_t i, size_t j) {
    const bool ni = std::isnan(a[i]), nj = std::isnan(a[j]);
    if (ni || nj) return !ni && nj;
    return a[i] < a[j];
  });
  size_t m = n;
  while (m > 0 && std::isnan(a[o[m - 1]])) out[o[--m]] = NAN;
  size_t i = 0;
  while (i < m) {
    size_t j = i;
    while (j + 1 < m && a[o[j + 1]] == a[o[i]]) ++j;
    const double r = (static_cast<double>(i) + static_cast<double>(j) + 2.0) / 2.0;
    for (size_t k = i; k <= j; ++k) out[o[k]] = r;
    i = j + 1;
  }
}

double spear(const double *x, const double *y, size_t n) {
  if (n < 2) return NAN;
  vector<double> rx(n), ry(n);
  midrank(x, n, rx.data());
  midrank(y, n, ry.data());
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i) {
    sx += rx[i];
    sy += ry[i];
    sxx += rx[i] * rx[i];
    syy += ry[i] * ry[i];
    sxy += rx[i] * ry[i];
  }
  const double dn = static_cast<double>(n), num = dn * sxy - sx * sy;
  const double den = (dn * sxx - sx * sx) * (dn * syy - sy * sy);
  if (den <= 0.0) return NAN;
  return num / std::sqrt(den);
}

double pears(const double *x, const double *y, size_t n) {
  if (n < 2) return NAN;
  const double mx = rmean(x, n), my = rmean(y, n);
  if (!std::isfinite(mx) || !std::isfinite(my)) return NAN;
  ld sxy = 0, sxx = 0, syy = 0;
  for (size_t i = 0; i < n; ++i) {
    const ld a = x[i] - mx, b = y[i] - my;
    sxy += a * b;
    sxx += a * a;
    syy += b * b;
  }
  if (sxx <= 0 || syy <= 0) return NAN;
  ld r = sxy / std::sqrt(sxx * syy);
  if (r > 1) r = 1;
  if (r < -1) r = -1;
  return static_cast<double>(r);
}

// The eigenvalues (descending) and vectors (columns) of the symmetric p x p
// a, by cyclic Jacobi rotations.
void jacobi(vector<double> a, int p, vector<double> &val, vector<double> &vec) {
  vec.assign(static_cast<size_t>(p) * p, 0.0);
  for (int i = 0; i < p; ++i) vec[static_cast<size_t>(i) * p + i] = 1;
  auto A = [&](int i, int j) -> double & { return a[static_cast<size_t>(j) * p + i]; };
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0;
    for (int i = 0; i < p; ++i) {
      for (int j = i + 1; j < p; ++j) off += A(i, j) * A(i, j);
    }
    if (off == 0) break;
    for (int i = 0; i < p; ++i) {
      for (int j = i + 1; j < p; ++j) {
        if (A(i, j) == 0) continue;
        const double th = (A(j, j) - A(i, i)) / (2 * A(i, j));
        const double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < p; ++k) {
          const double aki = A(k, i), akj = A(k, j);
          A(k, i) = c * aki - s * akj;
          A(k, j) = s * aki + c * akj;
        }
        for (int k = 0; k < p; ++k) {
          const double aik = A(i, k), ajk = A(j, k);
          A(i, k) = c * aik - s * ajk;
          A(j, k) = s * aik + c * ajk;
        }
        for (int k = 0; k < p; ++k) {
          double &vi = vec[static_cast<size_t>(i) * p + k], &vj = vec[static_cast<size_t>(j) * p + k];
          const double a1 = vi, a2 = vj;
          vi = c * a1 - s * a2;
          vj = s * a1 + c * a2;
        }
      }
    }
  }
  vector<int> o(static_cast<size_t>(p));
  for (int i = 0; i < p; ++i) o[static_cast<size_t>(i)] = i;
  std::sort(o.begin(), o.end(), [&](int i, int j) { return A(i, i) > A(j, j); });
  val.resize(static_cast<size_t>(p));
  vector<double> v2(static_cast<size_t>(p) * p);
  for (int k = 0; k < p; ++k) {
    val[static_cast<size_t>(k)] = A(o[static_cast<size_t>(k)], o[static_cast<size_t>(k)]);
    for (int i = 0; i < p; ++i) {
      v2[static_cast<size_t>(k) * p + i] = vec[static_cast<size_t>(o[static_cast<size_t>(k)]) * p + i];
    }
  }
  vec = v2;
}

// The minimum eigenvalue is at most 1e-8 of the largest in magnitude.
bool lowrank(const vector<double> &s, int p) {
  vector<double> v, e;
  jacobi(s, p, v, e);
  double mx = 2.220446049250313e-16;
  for (double x : v) mx = std::max(mx, std::fabs(x));
  return v.empty() || v.back() <= 1e-8 * mx;
}

// .rmbl_make_pd(): symmetrised, its eigenvalues floored at ridge max(1,
// max |ev|).
void makepd(vector<double> &s, int p, double ridge) {
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j) {
      const double m = (s[static_cast<size_t>(j) * p + i] + s[static_cast<size_t>(i) * p + j]) / 2;
      s[static_cast<size_t>(j) * p + i] = s[static_cast<size_t>(i) * p + j] = m;
    }
  }
  vector<double> v, e;
  jacobi(s, p, v, e);
  double mx = 0;
  for (double x : v) mx = std::max(mx, std::fabs(x));
  const double fl = ridge * std::max(1.0, mx);
  if (v.back() >= fl) return;
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j < p; ++j) {
      double t = 0;
      for (int k = 0; k < p; ++k) {
        t += e[static_cast<size_t>(k) * p + i] * std::max(v[static_cast<size_t>(k)], fl) * e[static_cast<size_t>(k) * p + j];
      }
      s[static_cast<size_t>(j) * p + i] = t;
    }
  }
}

// solve(a): the inverse by LU with partial pivoting; false when singular
// (a zero pivot, or a reciprocal condition below the machine epsilon).
bool inv(const vector<double> &a0, int p, vector<double> &av) {
  vector<double> a = a0;
  vector<int> piv(static_cast<size_t>(p));
  auto A = [&](int i, int j) -> double & { return a[static_cast<size_t>(j) * p + i]; };
  double an = 0;
  for (int j = 0; j < p; ++j) {
    double s = 0;
    for (int i = 0; i < p; ++i) s += std::fabs(a0[static_cast<size_t>(j) * p + i]);
    an = std::max(an, s);
  }
  for (int k = 0; k < p; ++k) {
    int m = k;
    for (int i = k + 1; i < p; ++i) {
      if (std::fabs(A(i, k)) > std::fabs(A(m, k))) m = i;
    }
    piv[static_cast<size_t>(k)] = m;
    if (A(m, k) == 0) return false;
    if (m != k) {
      for (int j = 0; j < p; ++j) std::swap(A(k, j), A(m, j));
    }
    for (int i = k + 1; i < p; ++i) {
      A(i, k) /= A(k, k);
      for (int j = k + 1; j < p; ++j) A(i, j) -= A(i, k) * A(k, j);
    }
  }
  av.assign(static_cast<size_t>(p) * p, 0.0);
  for (int c = 0; c < p; ++c) {
    vector<double> b(static_cast<size_t>(p), 0.0);
    b[static_cast<size_t>(c)] = 1;
    for (int k = 0; k < p; ++k) std::swap(b[static_cast<size_t>(k)], b[static_cast<size_t>(piv[static_cast<size_t>(k)])]);
    for (int i = 0; i < p; ++i) {
      for (int k = 0; k < i; ++k) b[static_cast<size_t>(i)] -= A(i, k) * b[static_cast<size_t>(k)];
    }
    for (int i = p - 1; i >= 0; --i) {
      for (int k = i + 1; k < p; ++k) b[static_cast<size_t>(i)] -= A(i, k) * b[static_cast<size_t>(k)];
      b[static_cast<size_t>(i)] /= A(i, i);
    }
    for (int i = 0; i < p; ++i) av[static_cast<size_t>(c) * p + i] = b[static_cast<size_t>(i)];
  }
  double ni = 0;
  for (int j = 0; j < p; ++j) {
    double s = 0;
    for (int i = 0; i < p; ++i) s += std::fabs(av[static_cast<size_t>(j) * p + i]);
    ni = std::max(ni, s);
  }
  return 1 / (an * ni) >= 2.220446049250313e-16;
}

// The covariance (n - 1) of the columns of x (n x p) as the reference
// kernel computes it: double sums, the mean a double sum over n.
vector<double> covm(const vector<double> &x, size_t n, int p) {
  vector<double> mu(static_cast<size_t>(p)), out(static_cast<size_t>(p) * p, NAN);
  if (n < 2) return out;
  for (int j = 0; j < p; ++j) {
    double s = 0;
    for (size_t i = 0; i < n; ++i) s += x[static_cast<size_t>(j) * n + i];
    mu[static_cast<size_t>(j)] = s / static_cast<double>(n);
  }
  for (int j = 0; j < p; ++j) {
    for (int k = j; k < p; ++k) {
      double s = 0;
      for (size_t i = 0; i < n; ++i) {
        s += (x[static_cast<size_t>(j) * n + i] - mu[static_cast<size_t>(j)]) *
             (x[static_cast<size_t>(k) * n + i] - mu[static_cast<size_t>(k)]);
      }
      out[static_cast<size_t>(k) * p + j] = out[static_cast<size_t>(j) * p + k] = s / (static_cast<double>(n) - 1.0);
    }
  }
  return out;
}

// R's cov() of complete rows (ok all 1) or pairwise complete rows.
vector<double> rcov(const double *x, size_t n, int p, bool pairwise) {
  vector<double> out(static_cast<size_t>(p) * p, NAN);
  if (!pairwise) {
    vector<double> xm(static_cast<size_t>(p));
    for (int j = 0; j < p; ++j) xm[static_cast<size_t>(j)] = rmean(x + static_cast<size_t>(j) * n, n);
    for (int i = 0; i < p; ++i) {
      for (int j = 0; j <= i; ++j) {
        ld s = 0;
        for (size_t k = 0; k < n; ++k) {
          s += (x[static_cast<size_t>(i) * n + k] - xm[static_cast<size_t>(i)]) *
               (x[static_cast<size_t>(j) * n + k] - xm[static_cast<size_t>(j)]);
        }
        out[static_cast<size_t>(i) * p + j] = out[static_cast<size_t>(j) * p + i] = static_cast<double>(s / (n - 1));
      }
    }
    return out;
  }
  for (int i = 0; i < p; ++i) {
    for (int j = 0; j <= i; ++j) {
      const double *a = x + static_cast<size_t>(i) * n, *b = x + static_cast<size_t>(j) * n;
      ld xs = 0, ys = 0;
      long nob = 0;
      for (size_t k = 0; k < n; ++k) {
        if (std::isnan(a[k]) || std::isnan(b[k])) continue;
        ++nob;
        xs += a[k];
        ys += b[k];
      }
      double v = NAN;
      if (nob >= 2) {
        xs /= nob;
        ys /= nob;
        ld s = 0;
        for (size_t k = 0; k < n; ++k) {
          if (std::isnan(a[k]) || std::isnan(b[k])) continue;
          s += (a[k] - xs) * (b[k] - ys);
        }
        v = static_cast<double>(s / (nob - 1));
      }
      out[static_cast<size_t>(i) * p + j] = out[static_cast<size_t>(j) * p + i] = v;
    }
  }
  return out;
}

const double eps = 2.220446049250313e-16;

}  // namespace

extern "C" {

int clab_hbin(const double *x, int n, int bins, int *cnt) CLAB_TRY {
  if (bins < 1 || bins > kutil::maxcell || n < 0) return CLAB_EARG;
  for (int b = 0; b < bins; ++b) cnt[b] = 0;
  double lo = INFINITY, hi = -INFINITY;
  for (int i = 0; i < n; ++i) {
    lo = std::min(lo, x[i]);
    hi = std::max(hi, x[i]);
  }
  if (n == 0 || !(lo < hi)) return CLAB_EDOM;
  // seq(lo, hi, length.out = bins + 1) as seq.default() builds it
  const int nb = bins + 1;
  vector<double> br(static_cast<size_t>(nb));
  br[0] = lo;
  br[static_cast<size_t>(bins)] = hi;
  const double by = (hi - lo) / bins;
  for (int i = 1; i < bins; ++i) br[static_cast<size_t>(i)] = lo + i * by;
  // hist(): the breaks moved by 1e-7 of the median width (of the range
  // when there are three breaks or fewer, of the smallest positive width
  // when four or five), the first one down
  vector<double> h(static_cast<size_t>(bins));
  for (int i = 0; i < bins; ++i) h[static_cast<size_t>(i)] = br[static_cast<size_t>(i) + 1] - br[static_cast<size_t>(i)];
  double dd;
  if (nb > 5) {
    vector<double> s = h;
    std::sort(s.begin(), s.end());
    const size_t m = s.size();
    dd = m % 2 ? s[m / 2] : (s[m / 2 - 1] + s[m / 2]) / 2;
  } else if (nb <= 3) {
    dd = hi - lo;
  } else {
    dd = INFINITY;
    for (double v : h) {
      if (v > 0) dd = std::min(dd, v);
    }
  }
  const double di = 1e-7 * dd;
  vector<double> fb(br);
  fb[0] -= di;
  for (int i = 1; i < nb; ++i) fb[static_cast<size_t>(i)] += di;
  for (int i = 0; i < n; ++i) {
    const double v = x[i];
    if (!std::isfinite(v)) continue;
    int a = 0, b = nb - 1;
    if (fb[0] <= v && (v < fb[static_cast<size_t>(b)] || v == fb[static_cast<size_t>(b)])) {
      while (b - a >= 2) {
        const int m = (a + b) / 2;
        if (v > fb[static_cast<size_t>(m)]) a = m;
        else b = m;
      }
      cnt[a]++;
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_pcor(const double *x, int n, int p, int meth, int minprs, double *cor, int *np) CLAB_TRY {
  if (n < 0 || p < 2) return CLAB_EARG;
  int k = 0;
  vector<double> a, b;
  for (int i = 0; i < p; ++i) {
    for (int j = i + 1; j < p; ++j, ++k) {
      a.clear();
      b.clear();
      for (int r = 0; r < n; ++r) {
        const double u = x[static_cast<size_t>(i) * n + r], v = x[static_cast<size_t>(j) * n + r];
        if (std::isnan(u) || std::isnan(v)) continue;
        a.push_back(u);
        b.push_back(v);
      }
      np[k] = static_cast<int>(a.size());
      cor[k] = np[k] < minprs ? NAN : meth == 0 ? spear(a.data(), b.data(), a.size()) : pears(a.data(), b.data(), a.size());
    }
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_nrun(const int *na, int n, int minrun, int *st, int *en) CLAB_TRY {
  int k = 0;
  for (int i = 0; i < n;) {
    if (!na[i]) {
      ++i;
      continue;
    }
    int j = i;
    while (j + 1 < n && na[j + 1]) ++j;
    if (j - i + 1 >= minrun) {
      st[k] = i + 1;
      en[k] = j + 1;
      ++k;
    }
    i = j + 1;
  }
  return k;
} CLAB_CATCH(CLAB_ENOMEM)

void clab_msum(const int *na, int n, int p, double *out) CLAB_TRY {
  double tot = 0, any = 0, all = 0, comp = 0;
  for (int j = 0; j < p; ++j) {
    double c = 0;
    for (int i = 0; i < n; ++i) c += na[static_cast<size_t>(j) * n + i] != 0;
    tot += c;
    any += c > 0;
    all += c == n && n > 0;
  }
  for (int i = 0; i < n; ++i) {
    bool ok = true;
    for (int j = 0; j < p && ok; ++j) ok = !na[static_cast<size_t>(j) * n + i];
    comp += ok;
  }
  const double cells = std::max(1.0, static_cast<double>(n) * p);
  out[0] = n;
  out[1] = p;
  out[2] = tot;
  out[3] = 100 * tot / cells;
  out[4] = n == 0 ? 0 : comp;
  out[5] = n == 0 ? 0 : 100 * comp / n;
  out[6] = any;
  out[7] = all;
} CLAB_CATCHV

void clab_pstat(const double *x, int n, double *out) CLAB_TRY {
  for (int k = 0; k < 7; ++k) out[k] = NAN;
  if (n < 1) return;
  const size_t N = static_cast<size_t>(n);
  const double mean = rmean(x, N);
  out[0] = mean;
  if (n >= 2) {
    double s = 0;
    for (int i = 0; i < n; ++i) s = std::max(s, std::fabs(x[i] - mean));
    if (s == 0) {
      out[1] = 0;
    } else {
      ld z2 = 0, z3 = 0;
      for (int i = 0; i < n; ++i) {
        const ld z = (x[i] - mean) / s, zz = z * z;
        z2 += zz;
        z3 += zz * z;
      }
      out[1] = std::sqrt(static_cast<double>(z2 / (static_cast<ld>(n) - 1.0L)) * s * s);
      const ld m2 = z2 / n;
      if (n > 2 && z2 > 0) out[2] = static_cast<double>((z3 / n) / std::pow(m2, 1.5L));
    }
  }
  vector<double> v(x, x + N);
  std::sort(v.begin(), v.end());
  out[3] = q7s(v, 0.5);
  out[4] = mad(x, N);
  const double q1 = q7s(v, 0.25), q3 = q7s(v, 0.75), iqr = q3 - q1;
  out[5] = q1 - 1.5 * iqr;
  out[6] = q3 + 1.5 * iqr;
} CLAB_CATCHV

void clab_q7(const double *x, int n, const double *pr, int np, double *out) CLAB_TRY {
  vector<double> v(x, x + n);
  std::sort(v.begin(), v.end());
  for (int k = 0; k < np; ++k) out[k] = q7s(v, pr[k]);
} CLAB_CATCHV

double clab_lpchisq(double x, double df) CLAB_TRY {
  const double q = clab_pchisq(x, df, 0);
  if (q > 1e-280 || !(x > 0)) return std::log(q);
  // the continued fraction of the upper incomplete gamma ratio, on the log
  // scale: log Q(a, x / 2) = a log(x / 2) - x / 2 - lgamma(a) + log(cf)
  const double a = df / 2, z = x / 2, tiny = 1e-300;
  double b = z + 1 - a, c = 1 / tiny, d = 1 / b, h = d;
  for (int i = 1; i < 100000; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = an * d + b;
    if (std::fabs(d) < tiny) d = tiny;
    c = b + an / c;
    if (std::fabs(c) < tiny) c = tiny;
    d = 1 / d;
    const double del = d * c;
    h *= del;
    if (std::fabs(del - 1) < eps) break;
  }
  return a * std::log(z) - z - std::lgamma(a) + std::log(h);
} CLAB_CATCH(kutil::nan())

int clab_maha(const double *x, int m, int p, int robust, double *d2) CLAB_TRY {
  if (m < p + 1 || p < 1) return CLAB_EARG;
  const size_t M = static_cast<size_t>(m);
  vector<double> z(M * static_cast<size_t>(p));
  for (int j = 0; j < p; ++j) {
    const double *c = x + static_cast<size_t>(j) * M;
    double ctr, sc;
    if (robust) {
      ctr = med(c, M);
      sc = mad(c, M);
      if (!(sc > 0)) {
        // the classical scale, R's sd(), or 1
        const double mu = rmean(c, M);
        ld s = 0;
        for (size_t i = 0; i < M; ++i) {
          const double e = c[i] - mu;
          s += e * e;
        }
        const double sd = std::sqrt(static_cast<double>(s / (M - 1)));
        sc = sd > 0 ? sd : 1;
      }
    } else {
      ld s = 0;
      for (size_t i = 0; i < M; ++i) s += c[i];
      ctr = static_cast<double>(s / M);
      double lo = c[0], hi = c[0];
      for (size_t i = 0; i < M; ++i) {
        lo = std::min(lo, c[i]);
        hi = std::max(hi, c[i]);
      }
      sc = hi - lo;
      if (!(std::isfinite(sc) && sc > 0)) sc = 1;
    }
    for (size_t i = 0; i < M; ++i) z[static_cast<size_t>(j) * M + i] = (c[i] - ctr) / sc;
  }
  vector<double> s = covm(z, M, p), si;
  if (lowrank(s, p)) return CLAB_EDOM;
  makepd(s, p, 1e-8);
  if (!inv(s, p, si)) return CLAB_EDOM;
  for (size_t i = 0; i < M; ++i) {
    double t = 0;
    for (int b = 0; b < p; ++b) {
      double zs = 0;
      for (int a = 0; a < p; ++a) zs += z[static_cast<size_t>(a) * M + i] * si[static_cast<size_t>(b) * p + a];
      t += zs * z[static_cast<size_t>(b) * M + i];
    }
    d2[i] = t;
  }
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_mcar(const double *x, int n, int p, int maxit, double tol, double *out, double *mu, double *sig) CLAB_TRY {
  if (n < 2 || p < 2) return CLAB_EARG;
  const size_t N = static_cast<size_t>(n), P = static_cast<size_t>(p);
  auto X = [&](size_t i, size_t j) { return x[j * N + i]; };
  // the rank check: complete-case covariance, or pairwise
  vector<double> cc;
  for (size_t i = 0; i < N; ++i) {
    bool ok = true;
    for (size_t j = 0; j < P && ok; ++j) ok = !std::isnan(X(i, j));
    if (ok) {
      for (size_t j = 0; j < P; ++j) cc.push_back(X(i, j));
    }
  }
  const size_t nc = cc.size() / P;
  vector<double> s0;
  if (nc >= P + 1) {
    vector<double> cm(nc * P);
    for (size_t i = 0; i < nc; ++i) {
      for (size_t j = 0; j < P; ++j) cm[j * nc + i] = cc[i * P + j];
    }
    s0 = rcov(cm.data(), nc, p, false);
  } else {
    s0 = rcov(x, N, p, true);
  }
  for (double &v : s0) {
    if (!std::isfinite(v)) v = 0;
  }
  out[0] = 1;
  if (lowrank(s0, p)) return CLAB_EDOM;
  out[0] = 2;
  // EM: start at the pairwise moments
  vector<double> m(P), s = rcov(x, N, p, true);
  for (size_t j = 0; j < P; ++j) {
    vector<double> o;
    for (size_t i = 0; i < N; ++i) {
      if (!std::isnan(X(i, j))) o.push_back(X(i, j));
    }
    m[j] = rmean(o.data(), o.size());
    if (std::isnan(m[j])) return CLAB_EDOM;
  }
  for (double &v : s) {
    if (!std::isfinite(v)) v = 0;
  }
  for (size_t j = 0; j < P; ++j) {
    double &d = s[j * P + j];
    if (!(std::isfinite(d) && d > 0)) d = 1;
  }
  const double ridge = 1e-8;
  makepd(s, p, ridge);
  int it = 0;
  vector<double> t1(P), t2(P * P), xi(P), ci(P * P);
  for (int k = 1; k <= maxit; ++k) {
    if (kutil::stop(static_cast<double>(N) * P * P)) return CLAB_EINTR;
    it = k;
    std::fill(t1.begin(), t1.end(), 0.0);
    std::fill(t2.begin(), t2.end(), 0.0);
    for (size_t i = 0; i < N; ++i) {
      vector<int> o, mi;
      for (size_t j = 0; j < P; ++j) (std::isnan(X(i, j)) ? mi : o).push_back(static_cast<int>(j));
      std::fill(ci.begin(), ci.end(), 0.0);
      if (mi.empty()) {
        for (size_t j = 0; j < P; ++j) xi[j] = X(i, j);
      } else if (o.empty()) {
        xi = m;
        ci = s;
      } else {
        const int no = static_cast<int>(o.size()), nm = static_cast<int>(mi.size());
        vector<double> soo(static_cast<size_t>(no) * no), sinv;
        for (int a = 0; a < no; ++a) {
          for (int b = 0; b < no; ++b) {
            soo[static_cast<size_t>(b) * no + a] = s[static_cast<size_t>(o[b]) * P + o[a]];
          }
        }
        if (!inv(soo, no, sinv)) {
          double md = 0;
          for (int a = 0; a < no; ++a) md += soo[static_cast<size_t>(a) * no + a];
          md /= no;
          for (int a = 0; a < no; ++a) soo[static_cast<size_t>(a) * no + a] += ridge * std::max(1.0, md);
          if (!inv(soo, no, sinv)) return CLAB_EDOM;
        }
        // B = S[m, o] S[o, o]^-1
        vector<double> bm(static_cast<size_t>(nm) * no, 0.0);
        for (int a = 0; a < nm; ++a) {
          for (int b = 0; b < no; ++b) {
            double t = 0;
            for (int c = 0; c < no; ++c) t += s[static_cast<size_t>(o[c]) * P + mi[a]] * sinv[static_cast<size_t>(b) * no + c];
            bm[static_cast<size_t>(b) * nm + a] = t;
          }
        }
        for (int b = 0; b < no; ++b) xi[o[b]] = X(i, o[b]);
        for (int a = 0; a < nm; ++a) {
          double t = 0;
          for (int b = 0; b < no; ++b) t += bm[static_cast<size_t>(b) * nm + a] * (X(i, o[b]) - m[o[b]]);
          xi[mi[a]] = m[mi[a]] + t;
        }
        for (int a = 0; a < nm; ++a) {
          for (int c = 0; c < nm; ++c) {
            double t = 0;
            for (int b = 0; b < no; ++b) t += bm[static_cast<size_t>(b) * nm + a] * s[static_cast<size_t>(mi[c]) * P + o[b]];
            ci[static_cast<size_t>(mi[c]) * P + mi[a]] = s[static_cast<size_t>(mi[c]) * P + mi[a]] - t;
          }
        }
      }
      for (size_t a = 0; a < P; ++a) {
        t1[a] += xi[a];
        for (size_t b = 0; b < P; ++b) t2[b * P + a] = (t2[b * P + a] + xi[a] * xi[b]) + ci[b * P + a];
      }
    }
    vector<double> mn(P), sn(P * P);
    for (size_t a = 0; a < P; ++a) mn[a] = t1[a] / n;
    for (size_t a = 0; a < P; ++a) {
      for (size_t b = 0; b < P; ++b) sn[b * P + a] = t2[b * P + a] / n - mn[a] * mn[b];
    }
    makepd(sn, p, ridge);
    double del = 0;
    for (size_t a = 0; a < P; ++a) del = std::max(del, std::fabs(mn[a] - m[a]));
    for (size_t a = 0; a < P * P; ++a) del = std::max(del, std::fabs(sn[a] - s[a]));
    m = mn;
    s = sn;
    if (std::isfinite(del) && del < tol) break;
  }
  // Little's d2 over the patterns, in the order they first appear
  std::map<std::string, int> seen;
  vector<std::string> pat;
  vector<vector<size_t>> rows;
  for (size_t i = 0; i < N; ++i) {
    std::string c(P, '0');
    for (size_t j = 0; j < P; ++j) c[j] = std::isnan(X(i, j)) ? '0' : '1';
    auto f = seen.find(c);
    if (f == seen.end()) {
      seen[c] = static_cast<int>(pat.size());
      pat.push_back(c);
      rows.push_back({i});
    } else {
      rows[static_cast<size_t>(f->second)].push_back(i);
    }
  }
  double d2 = 0, dfs = 0, used = 0, skip = 0;
  for (size_t k = 0; k < pat.size(); ++k) {
    vector<int> o;
    for (size_t j = 0; j < P; ++j) {
      if (pat[k][j] == '1') o.push_back(static_cast<int>(j));
    }
    const int no = static_cast<int>(o.size());
    const size_t nj = rows[k].size();
    vector<double> df(static_cast<size_t>(no)), sj(static_cast<size_t>(no) * no), si;
    for (int a = 0; a < no; ++a) {
      ld t = 0;
      for (size_t r : rows[k]) t += X(r, static_cast<size_t>(o[a]));
      df[static_cast<size_t>(a)] = static_cast<double>(t / nj) - m[static_cast<size_t>(o[a])];
      for (int b = 0; b < no; ++b) sj[static_cast<size_t>(b) * no + a] = s[static_cast<size_t>(o[b]) * P + o[a]];
    }
    if (!inv(sj, no, si)) {
      skip += 1;
      continue;
    }
    double q = 0;
    for (int b = 0; b < no; ++b) {
      double t = 0;
      for (int a = 0; a < no; ++a) t += df[static_cast<size_t>(a)] * si[static_cast<size_t>(b) * no + a];
      q += t * df[static_cast<size_t>(b)];
    }
    d2 += static_cast<double>(nj) * q;
    dfs += no;
    used += static_cast<double>(nj);
  }
  const double dof = dfs - p;
  out[0] = d2;
  out[1] = dof;
  out[2] = dof > 0 ? clab_pchisq(d2, dof, 0) : NAN;
  out[3] = dof > 0 ? clab_lpchisq(d2, dof) : NAN;
  out[4] = static_cast<double>(pat.size());
  out[5] = used;
  out[6] = it;
  out[7] = skip;
  for (size_t a = 0; a < P; ++a) mu[a] = m[a];
  for (size_t a = 0; a < P * P; ++a) sig[a] = s[a];
  return CLAB_OK;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_freq(const int *code, int n, int nlev, int sort, int maxlvl, int fold, int *lev, int *cnt) CLAB_TRY {
  vector<int> c(static_cast<size_t>(nlev), 0), o(static_cast<size_t>(nlev));
  for (int i = 0; i < n; ++i) {
    if (code[i] >= 1 && code[i] <= nlev) c[static_cast<size_t>(code[i] - 1)]++;
  }
  for (int i = 0; i < nlev; ++i) o[static_cast<size_t>(i)] = i;
  if (sort) std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return c[static_cast<size_t>(a)] > c[static_cast<size_t>(b)]; });
  int k = 0, other = 0;
  for (int i = 0; i < nlev; ++i) {
    if (i < maxlvl) {
      lev[k] = o[static_cast<size_t>(i)] + 1;
      cnt[k++] = c[static_cast<size_t>(o[static_cast<size_t>(i)])];
    } else {
      other += c[static_cast<size_t>(o[static_cast<size_t>(i)])];
    }
  }
  if (fold && nlev > maxlvl) {
    lev[k] = nlev + 1;
    cnt[k++] = other;
  }
  return k;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_mgrid(const int *na, int n, int p, int bins, int *start, int *shade) CLAB_TRY {
  // unique(as.integer(round(seq(0, n, length.out = bins + 1))))
  vector<int> e;
  const double by = static_cast<double>(n) / bins;
  for (int i = 0; i <= bins; ++i) {
    const double v = i == 0 ? 0.0 : i == bins ? static_cast<double>(n) : i * by;
    const int r = static_cast<int>(std::nearbyint(v));
    if (std::find(e.begin(), e.end(), r) == e.end()) e.push_back(r);
  }
  const int nb = static_cast<int>(e.size()) - 1;
  for (int b = 0; b < nb; ++b) {
    const int lo = e[static_cast<size_t>(b)], hi = e[static_cast<size_t>(b) + 1];
    start[b] = lo + 1;
    for (int j = 0; j < p; ++j) {
      int m = 0;
      for (int i = lo; i < hi; ++i) m += na[static_cast<size_t>(j) * n + i] != 0;
      const double fr = static_cast<double>(m) / (hi - lo);
      shade[static_cast<size_t>(j) * nb + b] = static_cast<int>(std::nearbyint(fr * 4));
    }
  }
  return nb;
} CLAB_CATCH(CLAB_ENOMEM)

int clab_dups(const int *code, int n, int nlev, int *rows, int *cnt) CLAB_TRY {
  vector<int> c(static_cast<size_t>(nlev) + 1, 0);
  for (int i = 0; i < n; ++i) c[static_cast<size_t>(code[i])]++;
  vector<int> o;
  for (int i = 0; i < n; ++i) {
    if (c[static_cast<size_t>(code[i])] > 1) o.push_back(i);
  }
  std::stable_sort(o.begin(), o.end(), [&](int a, int b) { return code[a] < code[b]; });
  for (size_t k = 0; k < o.size(); ++k) {
    rows[k] = o[k] + 1;
    cnt[k] = c[static_cast<size_t>(code[o[k]])];
  }
  return static_cast<int>(o.size());
} CLAB_CATCH(CLAB_ENOMEM)

void clab_cochq(const double *y, const double *s, int g, double *out) CLAB_TRY {
  ld sw = 0, swy = 0;
  vector<double> w(static_cast<size_t>(g));
  for (int i = 0; i < g; ++i) {
    w[static_cast<size_t>(i)] = 1 / (s[i] * s[i]);
    sw += w[static_cast<size_t>(i)];
    swy += w[static_cast<size_t>(i)] * y[i];
  }
  const double pe = static_cast<double>(swy) / static_cast<double>(sw);
  ld q = 0;
  for (int i = 0; i < g; ++i) {
    const double e = y[i] - pe;
    q += w[static_cast<size_t>(i)] * (e * e);
  }
  const double qs = static_cast<double>(q), df = g - 1;
  out[0] = pe;
  out[1] = std::sqrt(1 / static_cast<double>(sw));
  out[2] = qs;
  out[3] = df;
  out[4] = clab_pchisq(qs, df, 0);
  out[5] = qs > 0 ? std::max(0.0, (qs - df) / qs) : 0;
} CLAB_CATCHV

}  // extern "C"
