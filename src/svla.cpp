// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, the shared part of the scaling and count-dynamics kernels:
// the node tree of svk.h, dense algebra, R's random numbers and the
// distribution functions, the argument readers and the kernel table.
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <limits>

#include "crimlab.h"
#include "kutil.h"
#include "svk.h"

namespace svk {

const Node *Node::at(const char *name) const {
  for (std::size_t i = 0; i < nm.size() && i < l.size(); ++i) {
    if (nm[i] == name) return l[i].k == NUL ? nullptr : &l[i];
  }
  return nullptr;
}

Node &Node::add(const char *name, Node v) {
  k = LIST;
  nm.emplace_back(name);
  l.push_back(std::move(v));
  return l.back();
}

Node num(double x) {
  Node n;
  n.k = Node::NUM;
  n.d.assign(1, x);
  return n;
}

Node ival(double x) {
  Node n = num(x);
  n.k = Node::INT;
  return n;
}

Node lgl(bool b) {
  Node n = num(b ? 1.0 : 0.0);
  n.k = Node::LGL;
  return n;
}

Node str(const std::string &x) {
  Node n;
  n.k = Node::STR;
  n.s.assign(1, x);
  return n;
}

Node vec(const std::vector<double> &x) {
  Node n;
  n.k = Node::NUM;
  n.d = x;
  return n;
}

Node ivec(const std::vector<double> &x) {
  Node n = vec(x);
  n.k = Node::INT;
  return n;
}

Node list() {
  Node n;
  n.k = Node::LIST;
  return n;
}

Node mnode(const Mat &m) {
  Node n = vec(m.v);
  n.dim = {m.r, m.c};
  return n;
}

Node anode(const std::vector<double> &x, int a, int b, int c) {
  Node n = vec(x);
  n.dim = {a, b, c};
  return n;
}

Mat mul(const Mat &a, const Mat &b) {
  Mat o(a.r, b.c);
  for (int j = 0; j < b.c; ++j) {
    for (int k = 0; k < a.c; ++k) {
      const double x = b(k, j);
      if (x == 0.0) continue;
      for (int i = 0; i < a.r; ++i) o(i, j) += a(i, k) * x;
    }
  }
  return o;
}

Mat tmul(const Mat &a, const Mat &b) {
  Mat o(a.c, b.c);
  for (int i = 0; i < a.c; ++i) {
    for (int j = 0; j < b.c; ++j) {
      double s = 0;
      for (int k = 0; k < a.r; ++k) s += a(k, i) * b(k, j);
      o(i, j) = s;
    }
  }
  return o;
}

Mat mult(const Mat &a, const Mat &b) { return mul(a, tr(b)); }

Mat tr(const Mat &a) {
  Mat o(a.c, a.r);
  for (int i = 0; i < a.r; ++i) {
    for (int j = 0; j < a.c; ++j) o(j, i) = a(i, j);
  }
  return o;
}

Mat eye(int n) {
  Mat o(n, n);
  for (int i = 0; i < n; ++i) o(i, i) = 1.0;
  return o;
}

namespace {

// The sign of column j of v (and of u, when given) such that the largest
// component of v's column is positive.
void sgn(Mat &v, int j, Mat *u) {
  int b = 0;
  for (int i = 1; i < v.r; ++i) {
    if (std::fabs(v(i, j)) > std::fabs(v(b, j)) * (1 + 1e-12)) b = i;
  }
  if (v(b, j) >= 0) return;
  for (int i = 0; i < v.r; ++i) v(i, j) = -v(i, j);
  if (u) {
    for (int i = 0; i < u->r; ++i) (*u)(i, j) = -(*u)(i, j);
  }
}

}  // namespace

void eigsym(const Mat &a0, std::vector<double> &val, Mat &vec) {
  const int n = a0.r;
  Mat a = a0;
  vec = eye(n);
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0, tot = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        tot += a(i, j) * a(i, j);
        if (i != j) off += a(i, j) * a(i, j);
      }
    }
    if (off <= 1e-32 * tot || off == 0.0) break;
    for (int p = 0; p < n - 1; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = a(p, q);
        if (apq == 0.0) continue;
        const double th = (a(q, q) - a(p, p)) / (2 * apq);
        const double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < n; ++k) {
          const double akp = a(k, p), akq = a(k, q);
          a(k, p) = c * akp - s * akq;
          a(k, q) = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = a(p, k), aqk = a(q, k);
          a(p, k) = c * apk - s * aqk;
          a(q, k) = s * apk + c * aqk;
        }
        for (int k = 0; k < n; ++k) {
          const double vkp = vec(k, p), vkq = vec(k, q);
          vec(k, p) = c * vkp - s * vkq;
          vec(k, q) = s * vkp + c * vkq;
        }
      }
    }
  }
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  std::stable_sort(o.begin(), o.end(), [&](int x, int y) { return a(x, x) > a(y, y); });
  Mat v2(n, n);
  val.assign(static_cast<std::size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    const int s = o[static_cast<std::size_t>(j)];
    val[static_cast<std::size_t>(j)] = a(s, s);
    for (int i = 0; i < n; ++i) v2(i, j) = vec(i, s);
    sgn(v2, j, nullptr);
  }
  vec = v2;
}

void svd(const Mat &a, Mat &u, std::vector<double> &d, Mat &v) {
  if (a.r < a.c) {
    svd(tr(a), v, d, u);
    for (int j = 0; j < static_cast<int>(d.size()); ++j) {
      sgn(v, j, &u);  // the sign rule is on the right vectors
    }
    return;
  }
  const int m = a.r, n = a.c;
  Mat w = a;
  Mat vv = eye(n);
  for (int sweep = 0; sweep < 100; ++sweep) {
    bool rot = false;
    for (int p = 0; p < n - 1; ++p) {
      for (int q = p + 1; q < n; ++q) {
        double al = 0, be = 0, ga = 0;
        for (int i = 0; i < m; ++i) {
          al += w(i, p) * w(i, p);
          be += w(i, q) * w(i, q);
          ga += w(i, p) * w(i, q);
        }
        if (ga == 0.0 || std::fabs(ga) <= 1e-15 * std::sqrt(al * be)) continue;
        rot = true;
        const double ze = (be - al) / (2 * ga);
        const double t = (ze >= 0 ? 1.0 : -1.0) / (std::fabs(ze) + std::sqrt(1 + ze * ze));
        const double c = 1 / std::sqrt(1 + t * t), s = c * t;
        for (int i = 0; i < m; ++i) {
          const double x = w(i, p), y = w(i, q);
          w(i, p) = c * x - s * y;
          w(i, q) = s * x + c * y;
        }
        for (int i = 0; i < n; ++i) {
          const double x = vv(i, p), y = vv(i, q);
          vv(i, p) = c * x - s * y;
          vv(i, q) = s * x + c * y;
        }
      }
    }
    if (!rot) break;
  }
  std::vector<double> nr(static_cast<std::size_t>(n));
  for (int j = 0; j < n; ++j) {
    double s = 0;
    for (int i = 0; i < m; ++i) s += w(i, j) * w(i, j);
    nr[static_cast<std::size_t>(j)] = std::sqrt(s);
  }
  std::vector<int> o(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) o[static_cast<std::size_t>(i)] = i;
  std::stable_sort(o.begin(), o.end(), [&](int x, int y) { return nr[static_cast<std::size_t>(x)] > nr[static_cast<std::size_t>(y)]; });
  u = Mat(m, n);
  v = Mat(n, n);
  d.assign(static_cast<std::size_t>(n), 0.0);
  for (int j = 0; j < n; ++j) {
    const int s = o[static_cast<std::size_t>(j)];
    const double dj = nr[static_cast<std::size_t>(s)];
    d[static_cast<std::size_t>(j)] = dj;
    for (int i = 0; i < n; ++i) v(i, j) = vv(i, s);
    for (int i = 0; i < m; ++i) u(i, j) = dj > 0 ? w(i, s) / dj : 0.0;
    sgn(v, j, &u);
  }
}

bool chol(const Mat &a, Mat &r) {
  const int n = a.r;
  r = Mat(n, n);
  for (int j = 0; j < n; ++j) {
    double s = a(j, j);
    for (int k = 0; k < j; ++k) s -= r(k, j) * r(k, j);
    if (!(s > 0)) return false;
    r(j, j) = std::sqrt(s);
    for (int i = j + 1; i < n; ++i) {
      double t = a(j, i);
      for (int k = 0; k < j; ++k) t -= r(k, j) * r(k, i);
      r(j, i) = t / r(j, j);
    }
  }
  return true;
}

bool solve(const Mat &a0, const Mat &b, Mat &x) {
  const int n = a0.r;
  Mat a = a0;
  x = b;
  std::vector<int> piv(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) {
      if (std::fabs(a(i, k)) > std::fabs(a(p, k))) p = i;
    }
    if (a(p, k) == 0.0) return false;
    if (p != k) {
      for (int j = 0; j < n; ++j) std::swap(a(k, j), a(p, j));
      for (int j = 0; j < x.c; ++j) std::swap(x(k, j), x(p, j));
    }
    for (int i = k + 1; i < n; ++i) {
      const double f = a(i, k) / a(k, k);
      if (f == 0.0) continue;
      for (int j = k; j < n; ++j) a(i, j) -= f * a(k, j);
      for (int j = 0; j < x.c; ++j) x(i, j) -= f * x(k, j);
    }
  }
  for (int j = 0; j < x.c; ++j) {
    for (int i = n - 1; i >= 0; --i) {
      double s = x(i, j);
      for (int k = i + 1; k < n; ++k) s -= a(i, k) * x(k, j);
      x(i, j) = s / a(i, i);
    }
  }
  double mx = 0, mn = std::numeric_limits<double>::infinity();
  for (int i = 0; i < n; ++i) {
    mx = std::max(mx, std::fabs(a(i, i)));
    mn = std::min(mn, std::fabs(a(i, i)));
  }
  return n == 0 || mn > mx * 1e-14;
}

bool inv(const Mat &a, Mat &x) { return solve(a, eye(a.r), x); }

bool lsq(const Mat &x0, const std::vector<double> &y0, std::vector<double> &b) {
  const int m = x0.r, n = x0.c;
  if (m < n) return false;
  Mat a = x0;
  std::vector<double> y = y0;
  double big = 0;
  for (double e : a.v) big = std::max(big, std::fabs(e));
  for (int k = 0; k < n; ++k) {
    double nr = 0;
    for (int i = k; i < m; ++i) nr += a(i, k) * a(i, k);
    nr = std::sqrt(nr);
    if (!(nr > 1e-7 * big)) return false;
    const double al = a(k, k) > 0 ? -nr : nr;
    std::vector<double> v(static_cast<std::size_t>(m - k));
    for (int i = k; i < m; ++i) v[static_cast<std::size_t>(i - k)] = a(i, k);
    v[0] -= al;
    double vn = 0;
    for (double e : v) vn += e * e;
    if (vn == 0.0) continue;
    for (int j = k; j < n; ++j) {
      double s = 0;
      for (int i = k; i < m; ++i) s += v[static_cast<std::size_t>(i - k)] * a(i, j);
      s = 2 * s / vn;
      for (int i = k; i < m; ++i) a(i, j) -= s * v[static_cast<std::size_t>(i - k)];
    }
    double s = 0;
    for (int i = k; i < m; ++i) s += v[static_cast<std::size_t>(i - k)] * y[static_cast<std::size_t>(i)];
    s = 2 * s / vn;
    for (int i = k; i < m; ++i) y[static_cast<std::size_t>(i)] -= s * v[static_cast<std::size_t>(i - k)];
  }
  b.assign(static_cast<std::size_t>(n), 0.0);
  for (int i = n - 1; i >= 0; --i) {
    double s = y[static_cast<std::size_t>(i)];
    for (int k = i + 1; k < n; ++k) s -= a(i, k) * b[static_cast<std::size_t>(k)];
    b[static_cast<std::size_t>(i)] = s / a(i, i);
  }
  return true;
}

Mat pinv(const Mat &a, double tol) {
  Mat u, v;
  std::vector<double> d;
  svd(a, u, d, v);
  Mat o(a.c, a.r);
  for (std::size_t k = 0; k < d.size(); ++k) {
    if (!(d[k] > tol)) continue;
    const int kk = static_cast<int>(k);
    for (int i = 0; i < a.c; ++i) {
      for (int j = 0; j < a.r; ++j) o(i, j) += v(i, kk) * u(j, kk) / d[k];
    }
  }
  return o;
}

double det(const Mat &a0) {
  const int n = a0.r;
  Mat a = a0;
  double s = 1;
  for (int k = 0; k < n; ++k) {
    int p = k;
    for (int i = k + 1; i < n; ++i) {
      if (std::fabs(a(i, k)) > std::fabs(a(p, k))) p = i;
    }
    if (a(p, k) == 0.0) return 0.0;
    if (p != k) {
      for (int j = 0; j < n; ++j) std::swap(a(k, j), a(p, j));
      s = -s;
    }
    s *= a(k, k);
    for (int i = k + 1; i < n; ++i) {
      const double f = a(i, k) / a(k, k);
      for (int j = k; j < n; ++j) a(i, j) -= f * a(k, j);
    }
  }
  return s;
}

// --- R's generators ------------------------------------------------------

rng::rng(int seed) {
  std::uint32_t s = static_cast<std::uint32_t>(seed);
  for (int j = 0; j < 50; ++j) s = 69069u * s + 1u;
  s = 69069u * s + 1u;  // the first word is the position, reset to 624
  for (int j = 0; j < 624; ++j) {
    s = 69069u * s + 1u;
    mt[j] = s;
  }
  mti = 624;
}

double rng::unif() {
  static const std::uint32_t mag[2] = {0x0u, 0x9908b0dfu};
  std::uint32_t y;
  if (mti >= 624) {
    int kk;
    for (kk = 0; kk < 624 - 397; ++kk) {
      y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
      mt[kk] = mt[kk + 397] ^ (y >> 1) ^ mag[y & 0x1u];
    }
    for (; kk < 623; ++kk) {
      y = (mt[kk] & 0x80000000u) | (mt[kk + 1] & 0x7fffffffu);
      mt[kk] = mt[kk + (397 - 624)] ^ (y >> 1) ^ mag[y & 0x1u];
    }
    y = (mt[623] & 0x80000000u) | (mt[0] & 0x7fffffffu);
    mt[623] = mt[396] ^ (y >> 1) ^ mag[y & 0x1u];
    mti = 0;
  }
  y = mt[mti++];
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680u;
  y ^= (y << 15) & 0xefc60000u;
  y ^= (y >> 18);
  const double x = static_cast<double>(y) * 2.3283064365386963e-10;
  const double i2 = 2.328306437080797e-10;
  if (x <= 0.0) return 0.5 * i2;
  if (1.0 - x <= 0.0) return 1.0 - 0.5 * i2;
  return x;
}

double rng::norm() {
  const double big = 134217728;  // 2^27
  double u = unif();
  u = static_cast<int>(big * u) + unif();
  return qnorm(u / big);
}

double rng::expo() {
  static double q[16];
  static bool done = false;
  if (!done) {
    // q[k] = sum_{i = 1}^{k + 1} log(2)^i / i!
    double t = 1, s = 0;
    for (int k = 0; k < 16; ++k) {
      t *= 0.693147180559945309417232121458 / (k + 1);
      s += t;
      q[k] = s;
    }
    q[15] = 1.0;
    done = true;
  }
  double a = 0;
  double u = unif();
  while (u <= 0.0 || u >= 1.0) u = unif();
  for (;;) {
    u += u;
    if (u > 1.0) break;
    a += q[0];
  }
  u -= 1.0;
  if (u <= q[0]) return a + u;
  int i = 0;
  double us = unif(), um = us;
  do {
    us = unif();
    if (um > us) um = us;
    ++i;
  } while (u > q[i]);
  return a + um * q[0];
}

double rng::gamma(double a, double scale) {
  const double sqrt32 = 5.656854, exp_m1 = 0.36787944117144233;
  const double q1 = 0.04166669, q2 = 0.02083148, q3 = 0.00801191, q4 = 0.00144121, q5 = -7.388e-5,
               q6 = 2.4511e-4, q7 = 2.424e-4;
  const double a1 = 0.3333333, a2 = -0.250003, a3 = 0.2000062, a4 = -0.1662921, a5 = 0.1423657,
               a6 = -0.1367177, a7 = 0.1233795;
  if (std::isnan(a) || std::isnan(scale)) return kutil::nan();
  if (a <= 0.0 || scale <= 0.0) {
    if (scale == 0.0 || a == 0.0) return 0.0;
    return kutil::nan();
  }
  if (!std::isfinite(a) || !std::isfinite(scale)) return std::numeric_limits<double>::infinity();
  double e, p, q, r, t, u, v, w, x, ret;
  if (a < 1.0) {
    e = 1.0 + exp_m1 * a;
    for (;;) {
      p = e * unif();
      if (p >= 1.0) {
        x = -std::log((e - p) / a);
        if (expo() >= (1.0 - a) * std::log(x)) break;
      } else {
        x = std::exp(std::log(p) / a);
        if (expo() >= x) break;
      }
    }
    return scale * x;
  }
  if (a != aa) {
    aa = a;
    s2 = a - 0.5;
    s = std::sqrt(s2);
    dd = sqrt32 - s * 12.0;
  }
  t = norm();
  x = s + 0.5 * t;
  ret = x * x;
  if (t >= 0.0) return scale * ret;
  u = unif();
  if (dd * u <= t * t * t) return scale * ret;
  if (a != aaa) {
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r + q1) * r;
    if (a <= 3.686) {
      b = 0.463 + s + 0.178 * s2;
      si = 1.235;
      c = 0.195 / s - 0.079 + 0.16 * s;
    } else if (a <= 13.022) {
      b = 1.654 + 0.0076 * s2;
      si = 1.68 / s + 0.275;
      c = 0.062 / s + 0.024;
    } else {
      b = 1.77;
      si = 0.75;
      c = 0.1515 / s;
    }
  }
  if (x > 0.0) {
    v = t / (s + s);
    if (std::fabs(v) <= 0.25) {
      q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
    } else {
      q = q0 - s * t + 0.25 * t * t + (s2 + s2) * std::log(1.0 + v);
    }
    if (std::log(1.0 - u) <= q) return scale * ret;
  }
  for (;;) {
    e = expo();
    u = unif();
    u = u + u - 1.0;
    t = u < 0.0 ? b - si * e : b + si * e;
    if (t >= -0.71874483771719) {
      v = t / (s + s);
      if (std::fabs(v) <= 0.25) {
        q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
      } else {
        q = q0 - s * t + 0.25 * t * t + (s2 + s2) * std::log(1.0 + v);
      }
      if (q > 0.0) {
        w = std::expm1(q);
        if (c * std::fabs(u) <= w * std::exp(e - 0.5 * t * t)) break;
      }
    }
  }
  x = s + 0.5 * t;
  return scale * x * x;
}

double rng::index(double dn) {
  if (dn <= 0) return 0.0;
  const int bits = static_cast<int>(std::ceil(std::log2(dn)));
  double dv;
  do {
    std::int64_t vv = 0;
    for (int n = 0; n <= bits; n += 16) {
      const int v1 = static_cast<int>(std::floor(unif() * 65536));
      vv = 65536 * vv + v1;
    }
    dv = static_cast<double>(vv & ((static_cast<std::int64_t>(1) << bits) - 1));
  } while (dn <= dv);
  return dv;
}

// --- distributions -------------------------------------------------------

double pnorm(double x, bool lg) {
  if (!lg) return clab_pnorm(x);
  if (x > 0) return std::log1p(-clab_pnorm(-x));
  if (x > -37) return std::log(clab_pnorm(x));
  // the Mills-ratio series of log Phi(x) far in the lower tail
  const double w = 1 / (x * x);
  return -0.5 * x * x - std::log(-x) - 0.91893853320467274178 +
         std::log(1 - w + 3 * w * w - 15 * w * w * w + 105 * w * w * w * w);
}

double qnorm(double p) { return clab_qnorm(p); }

double dnorm(double x, bool lg) {
  if (std::isnan(x)) return x;
  if (lg) return -(0.918938533204672741780329736406 + 0.5 * x * x);
  x = std::fabs(x);
  if (x >= 2 * std::sqrt(std::numeric_limits<double>::max())) return 0.0;
  if (x < 5) return 0.398942280401432677939946059934 * std::exp(-0.5 * x * x);
  if (x > std::sqrt(-2 * 0.693147180559945309417232121458 * (std::numeric_limits<double>::min_exponent + 1 - std::numeric_limits<double>::digits))) {
    return 0.0;
  }
  const double x1 = std::ldexp(std::nearbyint(std::ldexp(x, 16)), -16), x2 = x - x1;
  return 0.398942280401432677939946059934 * (std::exp(-0.5 * x1 * x1) * std::exp((-0.5 * x2 - x1) * x2));
}

double pgamma(double x, double shape, double rate) { return clab_pgamma(x * rate, shape, 1); }

double qgamma(double p, double shape, double rate) {
  if (!(p > 0)) return 0.0;
  if (!(p < 1)) return std::numeric_limits<double>::infinity();
  // bracket, then Newton steps on the log scale guarded by bisection
  double lo = 0, hi = std::max(1.0, shape);
  while (clab_pgamma(hi, shape, 1) < p) hi *= 2;
  double x = 0.5 * (lo + hi);
  for (int it = 0; it < 300; ++it) {
    const double f = clab_pgamma(x, shape, 1) - p;
    if (f == 0.0) break;
    if (f < 0) {
      lo = x;
    } else {
      hi = x;
    }
    const double dens = std::exp((shape - 1) * std::log(x) - x - std::lgamma(shape));
    double nx = dens > 0 ? x - f / dens : 0.5 * (lo + hi);
    if (!(nx > lo && nx < hi)) nx = 0.5 * (lo + hi);
    if (std::fabs(nx - x) <= 1e-16 * x) {
      x = nx;
      break;
    }
    x = nx;
  }
  return x / rate;
}

double quant(std::vector<double> x, double p) {
  const std::size_t n = x.size();
  if (n == 0) return kutil::nan();
  std::sort(x.begin(), x.end());
  const double index = 1 + static_cast<double>(n - 1) * p;
  const double lo = std::floor(index), hi = std::ceil(index);
  const double ql = x[static_cast<std::size_t>(lo) - 1], qh = x[static_cast<std::size_t>(hi) - 1];
  if (index > lo && qh != ql) {
    const double h = index - lo;
    return (1 - h) * ql + h * qh;
  }
  return ql;
}

double sdv(const double *x, int n, int step) {
  if (n < 2) return kutil::nan();
  double m = 0;
  for (int i = 0; i < n; ++i) m += x[static_cast<std::size_t>(i) * step];
  m /= n;
  double s = 0;
  for (int i = 0; i < n; ++i) {
    const double e = x[static_cast<std::size_t>(i) * step] - m;
    s += e * e;
  }
  return std::sqrt(s / (n - 1));
}

// --- arguments -----------------------------------------------------------

std::string msg(const char *f, ...) {
  char b[512];
  std::va_list ap;
  va_start(ap, f);
  std::vsnprintf(b, sizeof b, f, ap);
  va_end(ap);
  return b;
}

bool has(const Node &a, const char *nm) { return a.at(nm) != nullptr; }

namespace {

bool isnum(const Node *n) { return n && (n->k == Node::NUM || n->k == Node::INT || n->k == Node::LGL); }

}  // namespace

bool getm(const Node &a, const char *nm, Mat &m, std::string &err, bool na, bool opt) {
  const Node *n = a.at(nm);
  if (!n && opt) {
    m = Mat();
    return true;
  }
  if (!isnum(n) || n->d.empty() || n->dim.size() > 2) {
    err = msg("`%s` must be a numeric matrix.", nm);
    return false;
  }
  if (n->dim.size() == 2) {
    m = Mat(n->dim[0], n->dim[1]);
  } else {
    m = Mat(static_cast<int>(n->d.size()), 1);
  }
  m.v = n->d;
  for (double x : m.v) {
    if (std::isinf(x) || (!na && std::isnan(x))) {
      err = na ? msg("`%s` must hold finite numbers or NA.", nm)
               : msg("`%s` must hold finite numbers with no missing values.", nm);
      return false;
    }
  }
  return true;
}

bool getv(const Node &a, const char *nm, std::vector<double> &v, std::string &err, bool na, bool opt) {
  const Node *n = a.at(nm);
  if (!n && opt) {
    v.clear();
    return true;
  }
  if (!isnum(n) || n->d.empty()) {
    err = msg("`%s` must be a numeric vector.", nm);
    return false;
  }
  v = n->d;
  for (double x : v) {
    if (std::isinf(x) || (!na && std::isnan(x))) {
      err = na ? msg("`%s` must hold finite numbers or NA.", nm)
               : msg("`%s` must hold finite numbers with no missing values.", nm);
      return false;
    }
  }
  return true;
}

bool getd(const Node &a, const char *nm, double &x, std::string &err, double lo, double hi) {
  const Node *n = a.at(nm);
  if (!isnum(n) || n->d.size() != 1 || !(n->d[0] >= lo && n->d[0] <= hi)) {
    // an infinite bound is no bound: it is left out of the message
    err = std::isinf(lo) && std::isinf(hi) ? msg("`%s` must be one number.", nm)
          : std::isinf(lo)                 ? msg("`%s` must be one number of at most %g.", nm, hi)
          : std::isinf(hi)                 ? msg("`%s` must be one number of at least %g.", nm, lo)
                                           : msg("`%s` must be one number from %g to %g.", nm, lo, hi);
    return false;
  }
  x = n->d[0];
  return true;
}

bool geti(const Node &a, const char *nm, int &x, std::string &err, int lo) {
  const Node *n = a.at(nm);
  if (!isnum(n) || n->d.size() != 1 || !(n->d[0] >= lo && n->d[0] <= 2147483647.0) ||
      n->d[0] != std::floor(n->d[0])) {
    err = msg("`%s` must be one whole number of at least %d.", nm, lo);
    return false;
  }
  x = static_cast<int>(n->d[0]);
  return true;
}

bool getb(const Node &a, const char *nm, bool &x, std::string &err) {
  const Node *n = a.at(nm);
  if (!n || n->k != Node::LGL || n->d.size() != 1 || std::isnan(n->d[0])) {
    err = msg("`%s` must be TRUE or FALSE.", nm);
    return false;
  }
  x = n->d[0] != 0;
  return true;
}

// --- the table -----------------------------------------------------------

const Fn table[] = {
    {"votagr", "x", k_votagr},
    {"amscal", "Z", k_amscal},
    {"bbscal", "X ndim=2 minscale=8", k_bbscal},
    {"optcls", "votes ndim=1 maxit=500 nstart=10 seed=42", k_optcls},
    {"dblctr", "D", k_dblctr},
    {"cmdsf", "D ndim=2", k_cmdsf},
    {"smacf", "D ndim=2 maxit=300 tol=1e-6 w? init?", k_smacf},
    {"nmdsf", "D ndim=2 maxit=300 tol=1e-6", k_nmdsf},
    {"mdsfit", "ev", k_mdsfit},
    {"unfstr", "xr xs D w?", k_unfstr},
    {"mlsmu", "D ndim=2 maxit=200 tol=1e-6 nstart=5", k_mlsmu},
    {"smunf", "D ndim=2 maxit=300 tol=1e-6", k_smunf},
    {"ipsrec", "xr xs?", k_ipsrec},
    {"nomutl", "x zy zn beta=15 w?", k_nomutl},
    {"nompr", "x zy zn beta=15 w?", k_nompr},
    {"nomll", "votes x zy zn beta=15 w?", k_nomll},
    {"procrs", "X target", k_procrs},
    {"irtlik", "votes x alpha beta", k_irtlik},
    {"irtpst", "chain std=TRUE", k_irtpst},
    {"ordocl", "Y ndim=2 maxit=500", k_ordocl},
    {"anchv", "y V ncat=5", k_anchv},
    {"indscl", "D ndim=2 maxit=300 tol=1e-6", k_indscl},
    {"nrmvec", "x y", k_nrmvec},
    {"cutlns", "normals cuts xlim=-1,1", k_cutlns},
    {"dwnom", "votes ndim=2 maxit=100", k_dwnom},
    {"nombt", "votes x normals cuts nboot=100 seed=42", k_nombt},
    {"irtem", "votes ndim=1 maxit=500 tol=1e-6 xprior=0,1 bprior=0,25 start? conv=cor", k_irtem},
    {"scboot", "Z method=am nboot=200 seed=42", k_scboot},
    {"wrdfsh", "dtm maxit=100 tol=1e-6", k_wrdfsh},
    {"bamscl", "Z nsamp=1000 burn=200 polarity=1 seed=42", k_bamscl},
    {"bmdsf", "D ndim=2 nsamp=1000 burn=200 sigma=1 seed=42", k_bmdsf},
    {"bunfld", "D ndim=2 nsamp=1000 burn=200 seed=42", k_bunfld},
    {"cjrirt", "votes ndim=1 nsamp=1000 burn=200 bvar=25 start? seed=0 chains=TRUE", k_cjrirt},
    {"anomin", "votes ndim=1 nsamp=500 burn=100 seed=42 thin=1 lop=0.025 minvotes=20 polarity=1 fixalpha=FALSE",
     k_anomin},
    {"ordirt", "Y ndim=1 nsamp=500 burn=100 seed=42 lprec=0", k_ordirt},
    {"dynirt",
     "votes period nsamp=500 burn=100 seed=42 thin=1 tau2=1 c0=-1 d0=-1 e0=0 E0=1 a0=0 A0=0.1 b0=0 B0=0.1 anchor?",
     k_dynirt},
    {"arimfc", "y h=12", k_arimfc},
    {"sarfc", "y h=12 order=1,1,1 seasonal=0,1,1,12", k_sarfc},
    {"catcor", "units", k_catcor},
    {"ousim", "y npath=100 days=365 dt=1 seed=42", k_ousim},
    {"fpdens", "y ngrid=64 nstep=200", k_fpdens},
    {"lvfit", "y year?", k_lvfit},
};
const int ntable = static_cast<int>(sizeof table / sizeof table[0]);

const Fn *find(const char *name) {
  for (const Fn &f : table) {
    if (std::strcmp(f.name, name) == 0) return &f;
  }
  return nullptr;
}

bool run(const char *name, const Node &a, Node &out, std::string &err) {
  const Fn *f = find(name);
  if (!f) {
    err = msg("no kernel %s.", name);
    return false;
  }
  out = list();
  return f->f(a, out, err);
}

}  // namespace svk
