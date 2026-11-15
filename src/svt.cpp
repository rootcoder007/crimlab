// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, scaling and ideal points: vote agreement, Aldrich-McKelvey
// and Basic Space (blackbox) scaling, Optimal Classification and its
// ordered form, classical, SMACOF and nonmetric MDS, unfolding, INDSCAL,
// the NOMINATE utility, likelihood, DW-NOMINATE fit and bootstrap, EM
// binary IRT, Wordfish, Procrustes rotation, anchoring vignettes and the
// geometry of cutting lines (Armstrong et al. 2021, Analyzing Spatial
// Models of Choice and Judgment, 2nd ed.). Starting values are drawn from
// R's generator seeded as set.seed(seed), so a fit equals the one of the
// same method run in R after set.seed().
#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "kutil.h"
#include "svk.h"

namespace svk {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

std::size_t z(int i) { return static_cast<std::size_t>(i); }

// Euclidean distances between the rows of x.
Mat pdist(const Mat &x) {
  Mat d(x.r, x.r);
  for (int i = 0; i < x.r; ++i) {
    for (int j = i + 1; j < x.r; ++j) {
      double s = 0;
      for (int k = 0; k < x.c; ++k) {
        const double e = x(i, k) - x(j, k);
        s += e * e;
      }
      d(i, j) = d(j, i) = std::sqrt(s);
    }
  }
  return d;
}

Mat rnormm(rng &g, int r, int c, double f = 1.0) {
  Mat m(r, c);
  for (double &e : m.v) e = g.norm() * f;
  return m;
}

double mean(const std::vector<double> &v) {
  double s = 0;
  for (double e : v) s += e;
  return s / static_cast<double>(v.size());
}

Mat dcen(const Mat &d) {
  const int n = d.r;
  Mat a(n, n), h(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      a(i, j) = d(i, j) * d(i, j);
      h(i, j) = (i == j ? 1.0 : 0.0) - 1.0 / n;
    }
  }
  Mat b = mul(mul(h, a), h);
  for (double &e : b.v) e *= -0.5;
  return b;
}

bool square(const Mat &d, const char *nm, std::string &err) {
  if (d.r != d.c || d.r < 2) {
    err = msg("`%s` must be a square matrix with two rows at least.", nm);
    return false;
  }
  return true;
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// seq(a, b, length.out = n)
std::vector<double> seqn(double a, double b, int n) {
  std::vector<double> v(z(n));
  if (n == 1) {
    v[0] = a;
    return v;
  }
  const double by = (b - a) / (n - 1);
  for (int i = 0; i < n; ++i) v[z(i)] = a + i * by;
  v[z(n - 1)] = b;
  return v;
}

// findInterval(x, v): the number of v[i] <= x, v sorted.
int fint(double x, const std::vector<double> &v) {
  int k = 0;
  for (double e : v) k += e <= x;
  return k;
}

// The coordinates of `ndim` dimensions for rows points (all the
// configurations of one fit), within kutil::maxcell before any is made.
bool dims(double rows, int nd, std::string &err) {
  return cap(rows * nd, kutil::maxcell, "`ndim` times the points scaled", err);
}

bool intr(std::string &err, double w) {
  if (!kutil::stop(w)) return false;
  err.clear();
  return true;
}

// Pseudo-inverse of a small symmetric matrix through its eigenvalues,
// those with |lambda| <= 1e-4 dropped (the DSYEV step of BLACKB).
Mat pinvsym(const Mat &a) {
  std::vector<double> ev;
  Mat v;
  eigsym(a, ev, v);
  Mat o(a.r, a.r);
  for (int k = 0; k < a.r; ++k) {
    if (!(std::fabs(ev[z(k)]) > 1e-4)) continue;
    for (int i = 0; i < a.r; ++i) {
      for (int j = 0; j < a.r; ++j) o(i, j) += v(i, k) * v(j, k) / ev[z(k)];
    }
  }
  return o;
}

}  // namespace

// --- vote agreement --------------------------------------------------------

bool k_votagr(const Node &a, Node &out, std::string &err) {
  Mat m;
  if (!getm(a, "x", m, err, true)) return false;
  const int n = m.r;
  Mat ag = eye(n);
  double s = 0, c = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int den = 0, same = 0;
      for (int k = 0; k < m.c; ++k) {
        if (std::isnan(m(i, k)) || std::isnan(m(j, k))) continue;
        ++den;
        same += m(i, k) == m(j, k);
      }
      ag(i, j) = ag(j, i) = den ? static_cast<double>(same) / den : NaN;
      if (den) {
        s += ag(i, j);
        ++c;
      }
    }
  }
  out.add("agree", mnode(ag));
  out.add("mean", num(c > 0 ? s / c : NaN));
  out.add("n", ival(n));
  out.add("m", ival(m.c));
  return true;
}

// --- Aldrich-McKelvey ------------------------------------------------------

namespace {

bool amstim(const Mat &Z, std::vector<double> &zh, std::string &err) {
  const int q = Z.c;
  if (q < 2) {
    err = "Aldrich-McKelvey scaling needs at least two stimuli.";
    return false;
  }
  Mat M(q, q);
  int used = 0;
  for (int i = 0; i < Z.r; ++i) {
    bool na = false;
    double mu = 0;
    for (int j = 0; j < q; ++j) {
      na = na || std::isnan(Z(i, j));
      mu += Z(i, j);
    }
    if (na) continue;
    mu /= q;
    double vr = 0;
    for (int j = 0; j < q; ++j) vr += (Z(i, j) - mu) * (Z(i, j) - mu);
    if (vr == 0.0) continue;
    // I - X (X'X)^-1 X' with X = [1, z]
    double s0 = q, s1 = 0, s2 = 0;
    for (int j = 0; j < q; ++j) {
      s1 += Z(i, j);
      s2 += Z(i, j) * Z(i, j);
    }
    Mat xtx(2, 2), ix;
    xtx(0, 0) = s0;
    xtx(0, 1) = xtx(1, 0) = s1;
    xtx(1, 1) = s2;
    inv(xtx, ix);
    for (int u = 0; u < q; ++u) {
      for (int v = 0; v < q; ++v) {
        const double zu = Z(i, u), zv = Z(i, v);
        const double p = ix(0, 0) + ix(0, 1) * zv + ix(1, 0) * zu + ix(1, 1) * zu * zv;
        M(u, v) += (u == v ? 1.0 : 0.0) - p;
      }
    }
    ++used;
  }
  if (used < 1) {
    err = "Aldrich-McKelvey scaling needs at least one respondent who places every stimulus and not all at "
          "the same point.";
    return false;
  }
  Mat C(q, q);
  for (int i = 0; i < q; ++i) {
    for (int j = 0; j < q; ++j) C(i, j) = (i == j ? 1.0 : 0.0) - 1.0 / q;
  }
  std::vector<double> ev;
  Mat vc;
  eigsym(mul(mul(C, M), C), ev, vc);
  int best = -1;
  for (int k = 0; k < q; ++k) {
    double cs = 0;
    for (int i = 0; i < q; ++i) cs += vc(i, k);
    if (std::fabs(cs) > 1e-8 * std::sqrt(static_cast<double>(q))) continue;  // the constant
    if (best < 0 || ev[z(k)] < ev[z(best)]) best = k;
  }
  if (best < 0) best = 0;
  zh.assign(z(q), 0.0);
  for (int i = 0; i < q; ++i) zh[z(i)] = vc(i, best);
  const double m = mean(zh), sd = sdv(zh.data(), q);
  for (double &e : zh) e = (e - m) / sd;
  if (zh[0] > 0) {
    for (double &e : zh) e = -e;
  }
  return true;
}

bool am(const Mat &Z, Node &out, std::string &err) {
  std::vector<double> zh;
  if (!amstim(Z, zh, err)) return false;
  const int n = Z.r;
  std::vector<double> al(z(n)), be(z(n));
  double sb = 0;
  for (int i = 0; i < n; ++i) {
    std::vector<double> y;
    Mat x(0, 2);
    std::vector<double> xs;
    for (int j = 0; j < Z.c; ++j) {
      if (std::isnan(Z(i, j))) continue;
      y.push_back(Z(i, j));
      xs.push_back(zh[z(j)]);
    }
    if (y.size() < 2) {
      al[z(i)] = 0;
      be[z(i)] = 1;
    } else {
      x = Mat(static_cast<int>(y.size()), 2);
      for (std::size_t k = 0; k < y.size(); ++k) {
        x(static_cast<int>(k), 0) = 1;
        x(static_cast<int>(k), 1) = xs[k];
      }
      std::vector<double> b;
      if (!lsq(x, y, b)) {
        err = "a respondent's placements could not be regressed on the stimulus positions.";
        return false;
      }
      al[z(i)] = b[0];
      be[z(i)] = std::fabs(b[1]) > 1e-10 ? b[1] : 1e-10;
    }
    sb += std::fabs(be[z(i)]);
  }
  std::vector<double> w(z(n));
  for (int i = 0; i < n; ++i) w[z(i)] = std::fabs(be[z(i)]) / sb * n;
  out.add("zhat", vec(zh));
  out.add("alpha", vec(al));
  out.add("beta", vec(be));
  out.add("weights", vec(w));
  return true;
}

}  // namespace

bool k_amscal(const Node &a, Node &out, std::string &err) {
  Mat Z;
  if (!getm(a, "Z", Z, err, true)) return false;
  return am(Z, out, err);
}

// --- Basic Space (blackbox) ------------------------------------------------

namespace {

using Obs = std::vector<char>;

std::vector<double> bbpol(const Mat &X, const Obs &ob) {
  const int np = X.r, ny = X.c;
  std::vector<int> fm(z(np), ny + 1);
  for (int i = 0; i < np; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (!ob[z(j) * z(np) + z(i)]) {
        fm[z(i)] = j + 1;
        break;
      }
    }
  }
  Mat R(ny, ny);
  for (int j = 0; j < ny; ++j) {
    for (int jj = 0; jj <= j; ++jj) {
      double n = 0, sa = 0, sb = 0, sab = 0, saa = 0, sbb = 0;
      for (int i = 0; i < np; ++i) {
        if (!(ob[z(j) * z(np) + z(i)] && jj + 1 < fm[z(i)])) continue;
        const double x = X(i, j), y = X(i, jj);
        n += 1;
        sa += x;
        sb += y;
        sab += x * y;
        saa += x * x;
        sbb += y * y;
      }
      const double AA = n * sab - sa * sb, BB = n * saa - sa * sa, CC = n * sbb - sb * sb;
      R(j, jj) = R(jj, j) = BB * CC <= 0 ? 0.0 : AA / std::sqrt(BB * CC);
    }
  }
  int best = 0;
  double bs = -1;
  for (int j = 0; j < ny; ++j) {
    double s = 0;
    for (int k = 0; k < ny; ++k) s += std::fabs(R(j, k));
    if (s > bs) {
      bs = s;
      best = j;
    }
  }
  std::vector<double> ll(z(ny));
  for (int k = 0; k < ny; ++k) ll[z(k)] = R(best, k) <= 0 ? -1.0 : 1.0;
  const int h = (ny - 1) / 2;
  for (int jk = 0; jk < ny; ++jk) {
    for (int j = 0; j < ny; ++j) {
      int neg = 0;
      for (int k = 0; k < ny; ++k) neg += R(j, k) * ll[z(k)] * ll[z(j)] < 0;
      if (neg > h) ll[z(j)] = -ll[z(j)];
    }
  }
  return ll;
}

struct Reg {
  Mat W, X;
  double sse = 0;
};

Reg bbreg(const Mat &XS, const Obs &ob, const Mat &P) {
  const int np = XS.r, ny = XS.c, k = P.c;
  Reg r;
  r.W = Mat(ny, k);
  r.X = Mat(np, ny, NaN);
  for (int j = 0; j < ny; ++j) {
    std::vector<int> o;
    for (int i = 0; i < np; ++i) {
      if (ob[z(j) * z(np) + z(i)]) o.push_back(i);
    }
    const int m = static_cast<int>(o.size());
    Mat Pp(m, k), y(m, 1);
    for (int t = 0; t < m; ++t) {
      for (int c = 0; c < k; ++c) Pp(t, c) = P(o[z(t)], c);
      y(t, 0) = XS(o[z(t)], j);
    }
    const Mat w = mul(pinvsym(tmul(Pp, Pp)), tmul(Pp, y));
    for (int c = 0; c < k; ++c) r.W(j, c) = w(c, 0);
    const Mat f = mul(Pp, w);
    for (int t = 0; t < m; ++t) {
      const double e = f(t, 0) - y(t, 0);
      r.X(o[z(t)], j) = e;
      r.sse += e * e;
    }
  }
  return r;
}

struct Reg2 {
  Mat P, X;
  double sse = 0;
};

Reg2 bbreg2(const Mat &XS, const Obs &ob, const Mat &W, const Mat &P0, int nf) {
  const int np = XS.r, ny = XS.c;
  Reg2 r;
  r.P = P0;
  r.X = Mat(np, ny, NaN);
  for (int i = 0; i < np; ++i) {
    std::vector<int> o;
    for (int j = 0; j < ny; ++j) {
      if (ob[z(j) * z(np) + z(i)]) o.push_back(j);
    }
    const int m = static_cast<int>(o.size());
    Mat A(m, nf), y(m, 1);
    for (int t = 0; t < m; ++t) {
      for (int c = 0; c < nf; ++c) A(t, c) = W(o[z(t)], c);
      y(t, 0) = XS(i, o[z(t)]) - W(o[z(t)], nf);
    }
    const Mat v = mul(pinvsym(tmul(A, A)), tmul(A, y));
    for (int c = 0; c < nf; ++c) r.P(i, c) = v(c, 0);
    const Mat f = mul(A, v);
    for (int t = 0; t < m; ++t) {
      const double e = f(t, 0) + W(o[z(t)], nf) - XS(i, o[z(t)]);
      r.X(i, o[z(t)]) = e;
      r.sse += e * e;
    }
  }
  return r;
}

struct BB {
  Mat psi, w;
  std::vector<double> c, d;
  double sse = 0, svsum = 0;
};

BB bbfit(const Mat &X0, int nf) {
  const int np = X0.r, ny = X0.c;
  Obs ob(X0.v.size());
  for (std::size_t i = 0; i < X0.v.size(); ++i) ob[i] = !std::isnan(X0.v[i]);
  std::vector<double> dc(z(ny));
  for (int j = 0; j < ny; ++j) {
    double s = 0, n = 0;
    for (int i = 0; i < np; ++i) {
      if (ob[z(j) * z(np) + z(i)]) {
        s += X0(i, j);
        n += 1;
      }
    }
    dc[z(j)] = s / n;
  }
  Mat X = X0, XS = X0;
  std::vector<double> ll = bbpol(X, ob);
  Mat PS(np, nf + 1);
  for (int jjj = 1; jjj <= nf; ++jjj) {
    std::vector<double> wxb(z(np));
    for (int i = 0; i < np; ++i) {
      double s = 0, n = 0;
      for (int j = 0; j < ny; ++j) {
        if (!ob[z(j) * z(np) + z(i)]) continue;
        s += (X(i, j) + (jjj == 1 ? -dc[z(j)] : 0.0)) * ll[z(j)];
        n += 1;
      }
      wxb[z(i)] = s / n;
    }
    const double mw = mean(wxb);
    double xxk = 0;
    for (double e : wxb) xxk += e * e;
    xxk -= np * mw * mw;
    Mat XT(np, 2);
    for (int i = 0; i < np; ++i) {
      XT(i, 0) = wxb[z(i)] - mw;
      XT(i, 1) = 1;
    }
    for (int mm = 0; mm < 4; ++mm) {
      const Reg r1 = bbreg(XS, ob, XT);
      const Reg2 r2 = bbreg2(XS, ob, r1.W, XT, 1);
      for (std::size_t e = 0; e < X.v.size(); ++e) {
        if (ob[e]) X.v[e] = r2.X.v[e];
      }
      double pxb = 0, pxs = 0;
      for (int i = 0; i < np; ++i) pxb += r2.P(i, 0);
      pxb /= np;
      for (int i = 0; i < np; ++i) pxs += r2.P(i, 0) * r2.P(i, 0);
      pxs -= np * pxb * pxb;
      const double xc = pxs > 0 ? std::sqrt(xxk / pxs) : 1.0;
      for (int i = 0; i < np; ++i) {
        XT(i, 0) = (r2.P(i, 0) - pxb) * xc;
        PS(i, jjj - 1) = XT(i, 0);
      }
    }
    if (jjj == nf) {
      for (int i = 0; i < np; ++i) PS(i, nf) = 1;
      X = X0;
    }
    XS = X;
    if (jjj < nf) ll = bbpol(X, ob);
  }
  Mat W;
  for (int nn = 0; nn < 5; ++nn) {
    const Reg r1 = bbreg(XS, ob, PS);
    W = r1.W;
    const Reg2 r2 = bbreg2(XS, ob, W, PS, nf);
    PS = r2.P;
    for (int c = 0; c < nf; ++c) {
      double m = 0;
      for (int i = 0; i < np; ++i) m += PS(i, c);
      m /= np;
      for (int i = 0; i < np; ++i) PS(i, c) -= m;
    }
    if (std::fabs(r1.sse - r2.sse) < 0.01) break;
  }
  Mat Pk(np, nf), Wk(ny, nf);
  for (int c = 0; c < nf; ++c) {
    for (int i = 0; i < np; ++i) Pk(i, c) = PS(i, c);
    for (int j = 0; j < ny; ++j) Wk(j, c) = W(j, c);
  }
  Mat u, v;
  std::vector<double> d;
  svd(mult(Pk, Wk), u, d, v);
  BB b;
  b.psi = Mat(np, nf);
  b.w = Mat(ny, nf);
  b.d.assign(d.begin(), d.begin() + nf);
  for (int c = 0; c < nf; ++c) {
    const double sq = std::sqrt(b.d[z(c)]);
    for (int i = 0; i < np; ++i) b.psi(i, c) = u(i, c) * sq;
    for (int j = 0; j < ny; ++j) b.w(j, c) = v(j, c) * sq;
  }
  b.c.assign(z(ny), 0.0);
  for (int j = 0; j < ny; ++j) b.c[z(j)] = W(j, nf);
  double so = 0, s2 = 0, no = 0;
  for (int i = 0; i < np; ++i) {
    for (int j = 0; j < ny; ++j) {
      if (!ob[z(j) * z(np) + z(i)]) continue;
      double f = b.c[z(j)];
      for (int c = 0; c < nf; ++c) f += b.psi(i, c) * b.w(j, c);
      const double e = f - X0(i, j);
      b.sse += e * e;
      so += X0(i, j);
      s2 += X0(i, j) * X0(i, j);
      no += 1;
    }
  }
  b.svsum = s2 - so * so / no;
  return b;
}

bool bb(const Mat &X, int ndim, int minscale, Node &out, std::string &err) {
  const int n = X.r, p = X.c;
  if (p < 2) {
    err = "Blackbox scaling needs at least two issues (columns).";
    return false;
  }
  std::vector<int> keep;
  const int need = std::max(1, std::min(minscale, p));
  for (int i = 0; i < n; ++i) {
    int c = 0;
    for (int j = 0; j < p; ++j) c += !std::isnan(X(i, j));
    if (c >= need) keep.push_back(i);
  }
  const int nk = static_cast<int>(keep.size());
  if (nk <= ndim) {
    err = "Too few respondents answer at least `minscale` issues.";
    return false;
  }
  Mat Xk(nk, p);
  for (int t = 0; t < nk; ++t) {
    for (int j = 0; j < p; ++j) Xk(t, j) = X(keep[z(t)], j);
  }
  for (int j = 0; j < p; ++j) {
    int c = 0;
    for (int t = 0; t < nk; ++t) c += !std::isnan(Xk(t, j));
    if (!c) {
      err = "Every issue needs at least one response among the scaled respondents.";
      return false;
    }
  }
  const int q = std::min(ndim, std::min(p - 1, nk - 1));
  const BB f = bbfit(Xk, q);
  Mat Psi(n, q, NaN);
  for (int t = 0; t < nk; ++t) {
    for (int c = 0; c < q; ++c) Psi(keep[z(t)], c) = f.psi(t, c);
  }
  std::vector<double> e2(f.d);
  for (double &e : e2) e *= e;
  out.add("ideal", mnode(Psi));
  out.add("weights", mnode(f.w));
  out.add("eigen", vec(e2));
  out.add("singular", vec(f.d));
  out.add("fit", num(f.svsum > 0 ? 1 - f.sse / f.svsum : 0.0));
  out.add("cmeans", vec(f.c));
  out.add("ndim", ival(q));
  return true;
}

}  // namespace

bool k_bbscal(const Node &a, Node &out, std::string &err) {
  Mat X;
  int nd, ms;
  if (!getm(a, "X", X, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "minscale", ms, err, 1)) return false;
  return bb(X, nd, ms, out, err);
}

// --- Optimal Classification ------------------------------------------------

namespace {

// The fewest errors of one vote split at a single point of the projections
// pv (projection, 1 for yea), either side voting yea; a split falls only
// between distinct projections. Adds them to te and returns the errors of
// voting with the majority.
int ocerr(std::vector<std::pair<double, int>> &pv, int &te) {
  std::sort(pv.begin(), pv.end());
  const int n = static_cast<int>(pv.size());
  int ny = 0;
  for (const auto &e : pv) ny += e.second;
  const int nul = std::min(ny, n - ny);
  int best = nul, below = 0;  // below: yeas at or below the split
  for (int k = 1; k < n; ++k) {
    below += pv[z(k - 1)].second;
    if (pv[z(k - 1)].first == pv[z(k)].first) continue;
    const int up = below + (n - k - (ny - below));  // yea above the split
    best = std::min(best, std::min(up, n - up));
  }
  te += best;
  return nul;
}

}  // namespace

bool k_optcls(const Node &a, Node &out, std::string &err) {
  Mat V;
  int nd, mi, ns, seed;
  if (!getm(a, "votes", V, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 1) ||
      !geti(a, "nstart", ns, err, 1) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const int n = V.r, m = V.c;
  if (!dims(static_cast<double>(n) + m, nd, err)) return false;
  rng g(seed);
  double best = std::numeric_limits<double>::infinity();
  for (int r = 0; r < ns; ++r) {
    Mat x = rnormm(g, n, nd);
    Mat nr(m, nd);
    for (int k = 0; k < mi; ++k) {
      if (intr(err, static_cast<double>(n) * m * nd)) return false;
      for (int j = 0; j < m; ++j) {
        std::vector<double> ym(z(nd)), nm(z(nd));
        int ny = 0, nn = 0;
        for (int i = 0; i < n; ++i) {
          if (V(i, j) == 1) {
            ++ny;
            for (int d = 0; d < nd; ++d) ym[z(d)] += x(i, d);
          } else if (V(i, j) == 0) {
            ++nn;
            for (int d = 0; d < nd; ++d) nm[z(d)] += x(i, d);
          }
        }
        if (!ny || !nn) continue;
        double l = 0;
        for (int d = 0; d < nd; ++d) {
          ym[z(d)] = ym[z(d)] / ny - nm[z(d)] / nn;
          l += ym[z(d)] * ym[z(d)];
        }
        l = std::sqrt(l);
        if (l > 0) {
          for (int d = 0; d < nd; ++d) nr(j, d) = ym[z(d)] / l;
        }
      }
      Mat xn(n, nd);
      for (int i = 0; i < n; ++i) {
        std::vector<double> dir(z(nd));
        int nv = 0;
        for (int j = 0; j < m; ++j) {
          if (std::isnan(V(i, j))) continue;
          ++nv;
          for (int d = 0; d < nd; ++d) dir[z(d)] += V(i, j) == 1 ? nr(j, d) : -nr(j, d);
        }
        if (!nv) continue;
        for (int d = 0; d < nd; ++d) xn(i, d) = dir[z(d)] / nv;
      }
      x = xn;
    }
    int te = 0, ne = 0;
    std::vector<std::pair<double, int>> pv;
    for (int j = 0; j < m; ++j) {
      pv.clear();
      for (int i = 0; i < n; ++i) {
        if (std::isnan(V(i, j))) continue;
        double pj = 0;
        for (int d = 0; d < nd; ++d) pj += x(i, d) * nr(j, d);
        pv.emplace_back(pj, V(i, j) == 1);
      }
      ne += ocerr(pv, te);
    }
    if (te < best) {
      best = te;
      out = list();
      out.add("ideal", mnode(x));
      out.add("normals", mnode(nr));
      out.add("pre", num(ne > 0 ? static_cast<double>(ne - te) / ne : 0.0));
      out.add("errors", ival(te));
      out.add("null", ival(ne));
      out.add("ndim", ival(nd));
    }
  }
  return true;
}

// --- MDS ---------------------------------------------------------------------

bool k_dblctr(const Node &a, Node &out, std::string &err) {
  Mat D;
  if (!getm(a, "D", D, err) || !square(D, "D", err)) return false;
  out.add("B", mnode(dcen(D)));
  return true;
}

bool k_cmdsf(const Node &a, Node &out, std::string &err) {
  Mat D;
  int nd;
  if (!getm(a, "D", D, err) || !square(D, "D", err) || !geti(a, "ndim", nd, err, 1)) return false;
  const int n = D.r;
  if (nd > n) {
    err = msg("`ndim` must be at most the number of objects, %d.", n);
    return false;
  }
  const Mat B = dcen(D);
  std::vector<double> ev;
  Mat vc;
  eigsym(B, ev, vc);
  Mat X(n, nd);
  double ps = 0, ta = 0;
  for (int k = 0; k < nd; ++k) {
    const double p = std::max(ev[z(k)], 0.0);
    ps += p;
    for (int i = 0; i < n; ++i) X(i, k) = vc(i, k) * std::sqrt(p);
  }
  for (double e : ev) ta += std::fabs(e);
  const Mat dm = pdist(X);
  double nu = 0, de = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (!(D(i, j) > 0)) continue;
      nu += (dm(i, j) - D(i, j)) * (dm(i, j) - D(i, j));
      de += D(i, j) * D(i, j);
    }
  }
  out.add("coords", mnode(X));
  out.add("eigen", vec(std::vector<double>(ev.begin(), ev.begin() + nd)));
  out.add("stress", num(de > 0 ? std::sqrt(nu / de) : 0.0));
  out.add("fit", num(ta > 0 ? ps / ta : 0.0));
  out.add("B", mnode(B));
  return true;
}

namespace {

// One Guttman transform of SMACOF: X <- V^+ B(X) X, with its stress.
Mat smstep(const Mat &X, const Mat &D, const Mat &W, const Mat &Vi, double &stress) {
  const int n = X.r;
  const Mat dx = pdist(X);
  Mat B(n, n);
  for (int i = 0; i < n; ++i) {
    double rs = 0;
    for (int j = 0; j < n; ++j) {
      if (i == j || !(dx(i, j) > 1e-12)) continue;
      B(i, j) = -W(i, j) * D(i, j) / dx(i, j);
      rs += B(i, j);
    }
    B(i, i) = -rs;
  }
  const Mat Xn = mul(mul(Vi, B), X);
  const Mat d2 = pdist(Xn);
  stress = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      const double r = D(i, j) - d2(i, j);
      stress += W(i, j) * r * r;
    }
  }
  stress *= 0.5;
  return Xn;
}

// V = diag(rowSums(W)) - W and its pseudo-inverse at the tolerance of
// Armadillo's pinv(): max(n, n) * max(sv) * epsilon.
Mat vinv(const Mat &W) {
  const int n = W.r;
  Mat V(n, n);
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j) s += W(i, j);
    for (int j = 0; j < n; ++j) V(i, j) = (i == j ? s : 0.0) - W(i, j);
  }
  Mat u, v;
  std::vector<double> d;
  svd(V, u, d, v);
  return pinv(V, n * d[0] * std::numeric_limits<double>::epsilon());
}

struct Sm {
  Mat X;
  double stress;
  int it;
  bool conv;
  bool stop = false;  // interrupted
};

Sm smacof(const Mat &D, Mat W, Mat X, int mi, double tol) {
  const int n = D.r;
  for (int i = 0; i < n; ++i) W(i, i) = 0;
  const Mat Vi = vinv(W);
  const Mat dx = pdist(X);
  double st = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) st += W(i, j) * (D(i, j) - dx(i, j)) * (D(i, j) - dx(i, j));
  }
  st /= 2;
  int it = 0;
  for (it = 1; it <= mi; ++it) {
    if (kutil::stop(static_cast<double>(n) * n * X.c)) return {X, st, it, false, true};
    double sn;
    const Mat Xn = smstep(X, D, W, Vi, sn);
    const bool done = std::fabs(st - sn) < tol;
    X = Xn;
    st = sn;
    if (done) break;
  }
  if (it > mi) it = mi;
  return {X, st, it, it < mi};
}

}  // namespace

bool k_smacf(const Node &a, Node &out, std::string &err) {
  Mat D, W, X;
  int nd, mi;
  double tol;
  if (!getm(a, "D", D, err) || !square(D, "D", err) || !geti(a, "ndim", nd, err, 1) ||
      !geti(a, "maxit", mi, err, 1) || !getd(a, "tol", tol, err, 0, HUGE_VAL) || !getm(a, "w", W, err, false, true) ||
      !getm(a, "init", X, err, false, true)) {
    return false;
  }
  const int n = D.r;
  if (W.v.empty()) {
    W = Mat(n, n, 1.0);
  } else if (W.r != n || W.c != n) {
    err = msg("`w` must be a %d by %d matrix.", n, n);
    return false;
  }
  if (X.v.empty()) {
    if (!dims(n, nd, err)) return false;
    rng g(42);
    X = rnormm(g, n, nd);
  } else if (X.r != n) {
    err = msg("`init` must have %d rows.", n);
    return false;
  }
  const Sm s = smacof(D, W, X, mi, tol);
  if (s.stop) return intr(err, 0);
  out.add("coords", mnode(s.X));
  out.add("stress", num(s.stress));
  out.add("iter", ival(s.it));
  out.add("conv", lgl(s.conv));
  return true;
}

namespace {

// Pool-adjacent-violators as the port it follows: one pass that steps back
// after each merge.
std::vector<double> pava(const std::vector<double> &y) {
  std::vector<double> val(y), w(y.size(), 1.0);
  std::vector<int> sz(y.size(), 1);
  std::size_t i = 0;
  while (i + 1 < val.size()) {
    if (val[i] > val[i + 1]) {
      const double nw = w[i] + w[i + 1];
      val[i] = (w[i] * val[i] + w[i + 1] * val[i + 1]) / nw;
      w[i] = nw;
      sz[i] += sz[i + 1];
      val.erase(val.begin() + static_cast<long>(i) + 1);
      w.erase(w.begin() + static_cast<long>(i) + 1);
      sz.erase(sz.begin() + static_cast<long>(i) + 1);
      if (i > 0) --i;
    } else {
      ++i;
    }
  }
  std::vector<double> o;
  for (std::size_t k = 0; k < val.size(); ++k) o.insert(o.end(), z(sz[k]), val[k]);
  return o;
}

}  // namespace

bool k_nmdsf(const Node &a, Node &out, std::string &err) {
  Mat D;
  int nd, mi;
  double tol;
  if (!getm(a, "D", D, err) || !square(D, "D", err) || !geti(a, "ndim", nd, err, 1) ||
      !geti(a, "maxit", mi, err, 1) || !getd(a, "tol", tol, err, 0, HUGE_VAL)) {
    return false;
  }
  const int n = D.r;
  if (!dims(n, nd, err)) return false;
  rng g(42);
  Mat X = rnormm(g, n, nd);
  std::vector<int> ui, uj;
  std::vector<double> d0;
  for (int j = 0; j < n; ++j) {
    for (int i = 0; i < j; ++i) {
      ui.push_back(i);
      uj.push_back(j);
      d0.push_back(D(i, j));
    }
  }
  const std::size_t np = d0.size();
  std::vector<std::size_t> ord(np);
  for (std::size_t k = 0; k < np; ++k) ord[k] = k;
  std::stable_sort(ord.begin(), ord.end(), [&](std::size_t p, std::size_t q) { return d0[p] < d0[q]; });
  double stress = 0;
  int it = 0;
  for (it = 1; it <= mi; ++it) {
    if (intr(err, static_cast<double>(n) * n * nd)) return false;
    const Mat dx = pdist(X);
    std::vector<double> dc(np), so(np);
    for (std::size_t k = 0; k < np; ++k) dc[k] = dx(ui[k], uj[k]);
    for (std::size_t k = 0; k < np; ++k) so[k] = dc[ord[k]];
    const std::vector<double> ho = pava(so);
    std::vector<double> dh(np);
    for (std::size_t k = 0; k < np; ++k) dh[ord[k]] = ho[k];
    Mat Dh(n, n);
    for (std::size_t k = 0; k < np; ++k) Dh(ui[k], uj[k]) = Dh(uj[k], ui[k]) = dh[k];
    const Sm s = smacof(Dh, Mat(n, n, 1.0), X, 1, 1e-6);
    double den = 0, nu = 0;
    for (std::size_t k = 0; k < np; ++k) {
      den += dc[k] * dc[k];
      nu += (dc[k] - dh[k]) * (dc[k] - dh[k]);
    }
    stress = den > 0 ? std::sqrt(nu / den) : 0.0;
    double mx = 0;
    for (std::size_t e = 0; e < X.v.size(); ++e) mx = std::max(mx, std::fabs(s.X.v[e] - X.v[e]));
    X = s.X;
    if (mx < tol) break;
  }
  if (it > mi) it = mi;
  out.add("coords", mnode(X));
  out.add("stress", num(stress));
  out.add("iter", ival(it));
  out.add("conv", lgl(it < mi));
  return true;
}

bool k_mdsfit(const Node &a, Node &out, std::string &err) {
  std::vector<double> ev;
  if (!getv(a, "ev", ev, err)) return false;
  double t = 0;
  for (double e : ev) t += std::fabs(e);
  std::vector<double> f, c;
  if (t != 0.0) {
    double s = 0;
    for (double e : ev) {
      f.push_back(std::max(e, 0.0) / t);
      s += f.back();
      c.push_back(s);
    }
  }
  out.add("fit", vec(f));
  out.add("cumfit", vec(c));
  out.add("ev", vec(ev));
  return true;
}

// --- unfolding ---------------------------------------------------------------

namespace {

double ustress(const Mat &xr, const Mat &xs, const Mat &D, const Mat *W) {
  double s = 0;
  for (int i = 0; i < xr.r; ++i) {
    for (int j = 0; j < xs.r; ++j) {
      if (std::isnan(D(i, j))) continue;
      double d = 0;
      for (int k = 0; k < xr.c; ++k) d += (xr(i, k) - xs(j, k)) * (xr(i, k) - xs(j, k));
      d = std::sqrt(d);
      s += (W ? (*W)(i, j) : 1.0) * (d - D(i, j)) * (d - D(i, j));
    }
  }
  return s;
}

void ccen(Mat &x) {
  for (int k = 0; k < x.c; ++k) {
    double m = 0;
    for (int i = 0; i < x.r; ++i) m += x(i, k);
    m /= x.r;
    for (int i = 0; i < x.r; ++i) x(i, k) -= m;
  }
}

}  // namespace

bool k_unfstr(const Node &a, Node &out, std::string &err) {
  Mat xr, xs, D, W;
  if (!getm(a, "xr", xr, err) || !getm(a, "xs", xs, err) || !getm(a, "D", D, err, true) ||
      !getm(a, "w", W, err, false, true)) {
    return false;
  }
  if (xr.c != xs.c || D.r != xr.r || D.c != xs.r || (!W.v.empty() && (W.r != D.r || W.c != D.c))) {
    err = "`D` (and `w`) must have a row per row of `xr` and a column per row of `xs`, and `xr` and `xs` the "
          "same columns.";
    return false;
  }
  out.add("stress", num(ustress(xr, xs, D, W.v.empty() ? nullptr : &W)));
  return true;
}

bool k_mlsmu(const Node &a, Node &out, std::string &err) {
  Mat D;
  int nd, mi, ns;
  double tol;
  if (!getm(a, "D", D, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 1) ||
      !getd(a, "tol", tol, err, 0, HUGE_VAL) || !geti(a, "nstart", ns, err, 1)) {
    return false;
  }
  const int nr = D.r, nc = D.c;
  if (!dims(static_cast<double>(nr) + nc, nd, err)) return false;
  rng g(42);
  double best = std::numeric_limits<double>::infinity();
  bool any = false;
  for (int r = 0; r < ns; ++r) {
    Mat Xr = rnormm(g, nr, nd), Xs = rnormm(g, nc, nd);
    ccen(Xr);
    ccen(Xs);
    Mat Dh(nr, nc);
    for (int i = 0; i < nr; ++i) {
      double s = 0, c = 0;
      for (int j = 0; j < nc; ++j) {
        if (!std::isnan(D(i, j))) {
          s += D(i, j);
          c += 1;
        }
      }
      for (int j = 0; j < nc; ++j) Dh(i, j) = D(i, j) - s / c;
    }
    double prev = std::numeric_limits<double>::infinity(), stress = prev;
    int it = 0;
    for (it = 1; it <= mi; ++it) {
      if (intr(err, static_cast<double>(nr) * nc * nd)) return false;
      Mat dm(nr, nc);
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          double s = 0;
          for (int k = 0; k < nd; ++k) s += (Xr(i, k) - Xs(j, k)) * (Xr(i, k) - Xs(j, k));
          dm(i, j) = std::max(std::sqrt(s), 1e-12);
        }
      }
      Mat gr(nr, nd);
      for (int i = 0; i < nr; ++i) {
        for (int j = 0; j < nc; ++j) {
          for (int k = 0; k < nd; ++k) {
            gr(i, k) += 2 * (dm(i, j) - Dh(i, j)) * (Xr(i, k) - Xs(j, k)) / dm(i, j);
          }
        }
      }
      for (double &e : gr.v) e /= nc;
      std::vector<double> ev;
      Mat vv;
      eigsym(tmul(gr, gr), ev, vv);
      const double gmr = 2 / (nc * std::max(ev[0], 1e-10));
      for (std::size_t e = 0; e < Xr.v.size(); ++e) Xr.v[e] -= gmr * gr.v[e];
      ccen(Xr);
      Mat gs(nc, nd);
      for (int j = 0; j < nc; ++j) {
        for (int i = 0; i < nr; ++i) {
          for (int k = 0; k < nd; ++k) {
            gs(j, k) += 2 * (dm(i, j) - Dh(i, j)) * (Xs(j, k) - Xr(i, k)) / dm(i, j);
          }
        }
      }
      for (double &e : gs.v) e /= nr;
      eigsym(tmul(gs, gs), ev, vv);
      const double gms = 2 / (nr * std::max(ev[0], 1e-10));
      for (std::size_t e = 0; e < Xs.v.size(); ++e) Xs.v[e] -= gms * gs.v[e];
      ccen(Xs);
      stress = ustress(Xr, Xs, D, nullptr);
      if (std::isfinite(prev) && std::fabs(prev - stress) / std::max(prev, 1e-12) < tol) break;
      prev = stress;
    }
    if (it > mi) it = mi;
    if (stress < best) {
      best = stress;
      any = true;
      out = list();
      out.add("xr", mnode(Xr));
      out.add("xs", mnode(Xs));
      out.add("stress", num(stress));
      out.add("iter", ival(it));
      out.add("conv", lgl(it < mi));
    }
  }
  if (!any) {
    err = "no start reached a finite stress.";
    return false;
  }
  return true;
}

bool k_smunf(const Node &a, Node &out, std::string &err) {
  Mat D;
  int nd, mi;
  double tol;
  if (!getm(a, "D", D, err) || !geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 1) ||
      !getd(a, "tol", tol, err, 0, HUGE_VAL)) {
    return false;
  }
  const int nr = D.r, nc = D.c, n = nr + nc;
  if (!dims(2.0 * n, nd, err)) return false;
  rng g(42);
  const Mat Xr = rnormm(g, nr, nd), Xs = rnormm(g, nc, nd);
  Mat Df(n, n), W(n, n), X(n, nd);
  for (int i = 0; i < nr; ++i) {
    for (int j = 0; j < nc; ++j) {
      Df(i, nr + j) = Df(nr + j, i) = D(i, j);
      W(i, nr + j) = W(nr + j, i) = 1;
    }
  }
  for (int k = 0; k < nd; ++k) {
    for (int i = 0; i < nr; ++i) X(i, k) = Xr(i, k);
    for (int j = 0; j < nc; ++j) X(nr + j, k) = Xs(j, k);
  }
  std::vector<double> vi(z(n));
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < n; ++j) s += W(i, j);
    vi[z(i)] = s > 1e-12 ? 1 / s : 0.0;
  }
  Mat dx = pdist(X);
  double st = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) st += W(i, j) * (Df(i, j) - dx(i, j)) * (Df(i, j) - dx(i, j));
  }
  st /= 2;
  int it = 0;
  for (it = 1; it <= mi; ++it) {
    if (intr(err, static_cast<double>(n) * n * nd)) return false;
    Mat B(n, n);
    for (int i = 0; i < n; ++i) {
      double rs = 0;
      for (int j = 0; j < n; ++j) {
        if (i != j && dx(i, j) > 1e-12) B(i, j) = -W(i, j) * Df(i, j) / dx(i, j);
        rs += B(i, j);
      }
      B(i, i) = -rs;
    }
    Mat Xn = mul(B, X);
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < nd; ++k) Xn(i, k) *= vi[z(i)];
    }
    dx = pdist(Xn);
    double sn = 0;
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) sn += W(i, j) * (Df(i, j) - dx(i, j)) * (Df(i, j) - dx(i, j));
    }
    sn /= 2;
    const bool done = std::fabs(st - sn) < tol;
    X = Xn;
    st = sn;
    if (done) break;
  }
  if (it > mi) it = mi;
  Mat xr(nr, nd), xs(nc, nd);
  for (int k = 0; k < nd; ++k) {
    for (int i = 0; i < nr; ++i) xr(i, k) = X(i, k);
    for (int j = 0; j < nc; ++j) xs(j, k) = X(nr + j, k);
  }
  out.add("xr", mnode(xr));
  out.add("xs", mnode(xs));
  out.add("stress", num(st));
  out.add("iter", ival(it));
  out.add("conv", lgl(it < mi));
  return true;
}

bool k_ipsrec(const Node &a, Node &out, std::string &err) {
  Mat xr, xs;
  if (!getm(a, "xr", xr, err) || !getm(a, "xs", xs, err, false, true)) return false;
  out.add("ideal", mnode(xr));
  return true;
}

// --- NOMINATE ----------------------------------------------------------------

namespace {

bool nomargs(const Node &a, Mat &x, Mat &zy, Mat &zn, double &beta, std::vector<double> &w, std::string &err) {
  if (!getm(a, "x", x, err) || !getm(a, "zy", zy, err) || !getm(a, "zn", zn, err) ||
      !getd(a, "beta", beta, err, -HUGE_VAL, HUGE_VAL) || !getv(a, "w", w, err, false, true)) {
    return false;
  }
  if (zy.c != x.c || zn.c != x.c || zn.r != zy.r) {
    err = "`zy` and `zn` must have the same rows, and the columns of `x`.";
    return false;
  }
  if (w.empty()) w.assign(z(x.c), 1.0);
  if (static_cast<int>(w.size()) != x.c) {
    err = msg("`w` must have one weight per dimension, %d.", x.c);
    return false;
  }
  return true;
}

void nomu(const Mat &x, const Mat &zy, const Mat &zn, double beta, const std::vector<double> &w, Mat &uy, Mat &un,
          Mat &p) {
  uy = Mat(x.r, zy.r);
  un = Mat(x.r, zy.r);
  p = Mat(x.r, zy.r);
  for (int i = 0; i < x.r; ++i) {
    for (int j = 0; j < zy.r; ++j) {
      double dy = 0, dn = 0;
      for (int k = 0; k < x.c; ++k) {
        dy += w[z(k)] * w[z(k)] * (x(i, k) - zy(j, k)) * (x(i, k) - zy(j, k));
        dn += w[z(k)] * w[z(k)] * (x(i, k) - zn(j, k)) * (x(i, k) - zn(j, k));
      }
      uy(i, j) = beta * std::exp(-0.5 * dy);
      un(i, j) = beta * std::exp(-0.5 * dn);
      p(i, j) = 1 / (1 + std::exp(-(uy(i, j) - un(i, j))));
    }
  }
}

}  // namespace

bool k_nomutl(const Node &a, Node &out, std::string &err) {
  Mat x, zy, zn, uy, un, p;
  double beta;
  std::vector<double> w;
  if (!nomargs(a, x, zy, zn, beta, w, err)) return false;
  nomu(x, zy, zn, beta, w, uy, un, p);
  Mat df = uy;
  for (std::size_t e = 0; e < df.v.size(); ++e) df.v[e] -= un.v[e];
  out.add("uyea", mnode(uy));
  out.add("unay", mnode(un));
  out.add("udiff", mnode(df));
  out.add("prob", mnode(p));
  return true;
}

bool k_nompr(const Node &a, Node &out, std::string &err) {
  std::vector<double> x, zy, zn, w;
  double beta;
  if (!getv(a, "x", x, err) || !getv(a, "zy", zy, err) || !getv(a, "zn", zn, err) ||
      !getd(a, "beta", beta, err, -HUGE_VAL, HUGE_VAL) || !getv(a, "w", w, err, false, true)) {
    return false;
  }
  if (w.empty()) w.assign(x.size(), 1.0);
  if (zy.size() != x.size() || zn.size() != x.size() || w.size() != x.size()) {
    err = "`x`, `zy`, `zn` and `w` must have the same length.";
    return false;
  }
  double dy = 0, dn = 0;
  for (std::size_t k = 0; k < x.size(); ++k) {
    dy += w[k] * w[k] * (x[k] - zy[k]) * (x[k] - zy[k]);
    dn += w[k] * w[k] * (x[k] - zn[k]) * (x[k] - zn[k]);
  }
  const double v = beta * (std::exp(-0.5 * dy) - std::exp(-0.5 * dn));
  out.add("prob", num(1 / (1 + std::exp(-v))));
  return true;
}

bool k_nomll(const Node &a, Node &out, std::string &err) {
  Mat V, x, zy, zn, uy, un, p;
  double beta;
  std::vector<double> w;
  if (!getm(a, "votes", V, err, true) || !nomargs(a, x, zy, zn, beta, w, err)) return false;
  if (V.r != x.r || V.c != zy.r) {
    err = "`votes` must have a row per row of `x` and a column per row of `zy`.";
    return false;
  }
  nomu(x, zy, zn, beta, w, uy, un, p);
  double ll = 0;
  int nc = 0, nt = 0;
  for (int i = 0; i < V.r; ++i) {
    for (int j = 0; j < V.c; ++j) {
      if (std::isnan(V(i, j))) continue;
      const double q = std::min(std::max(p(i, j), 1e-10), 1 - 1e-10);
      if (V(i, j) == 1) {
        ll += std::log(q);
        nc += q > 0.5;
      } else {
        ll += std::log(1 - q);
        nc += q < 0.5;
      }
      ++nt;
    }
  }
  out.add("loglik", num(ll));
  out.add("gmp", num(nt > 0 ? static_cast<double>(nc) / nt : 0.0));
  out.add("correct", ival(nc));
  out.add("total", ival(nt));
  return true;
}

// --- Procrustes --------------------------------------------------------------

bool k_procrs(const Node &a, Node &out, std::string &err) {
  Mat X, T;
  if (!getm(a, "X", X, err) || !getm(a, "target", T, err)) return false;
  if (X.r != T.r || X.c != T.c) {
    err = "`X` and `target` must have the same dimensions.";
    return false;
  }
  const int n = X.r, d = X.c;
  std::vector<double> mt(z(d));
  Mat Xc = X, Tc = T;
  ccen(Xc);
  for (int k = 0; k < d; ++k) {
    for (int i = 0; i < n; ++i) mt[z(k)] += T(i, k);
    mt[z(k)] /= n;
  }
  ccen(Tc);
  Mat u, v;
  std::vector<double> s;
  svd(tmul(Xc, Tc), u, s, v);
  Mat Q = mult(u, v);
  if (det(Q) < 0) {
    for (int i = 0; i < v.r; ++i) v(i, v.c - 1) = -v(i, v.c - 1);
    Q = mult(u, v);
  }
  Mat R = mul(Xc, Q);
  double mse = 0, sc = 0;
  for (int k = 0; k < d; ++k) {
    for (int i = 0; i < n; ++i) {
      R(i, k) += mt[z(k)];
      mse += (R(i, k) - T(i, k)) * (R(i, k) - T(i, k));
    }
  }
  for (double e : s) sc += e;
  out.add("rotated", mnode(R));
  out.add("rotation", mnode(Q));
  out.add("scale", num(sc));
  out.add("mse", num(mse / (static_cast<double>(n) * d)));
  return true;
}

// --- IRT likelihood and posterior summaries ---------------------------------

bool k_irtlik(const Node &a, Node &out, std::string &err) {
  Mat V, x, b;
  std::vector<double> al;
  if (!getm(a, "votes", V, err, true) || !getm(a, "x", x, err) || !getv(a, "alpha", al, err) ||
      !getm(a, "beta", b, err)) {
    return false;
  }
  if (x.r != V.r || b.r != V.c || b.c != x.c || static_cast<int>(al.size()) != V.c) {
    err = "`x` must have a row per legislator, `alpha` and `beta` one entry (row) per vote, and `beta` the "
          "columns of `x`.";
    return false;
  }
  Mat P(V.r, V.c);
  double ll = 0;
  int nc = 0, nt = 0;
  for (int i = 0; i < V.r; ++i) {
    for (int j = 0; j < V.c; ++j) {
      if (std::isnan(V(i, j))) continue;
      double zz = 0;
      for (int k = 0; k < x.c; ++k) zz += b(j, k) * x(i, k);
      zz -= al[z(j)];
      const double p = std::min(std::max(pnorm(zz), 1e-10), 1 - 1e-10);
      P(i, j) = p;
      if (V(i, j) == 1) {
        ll += std::log(p);
        nc += p > 0.5;
      } else {
        ll += std::log(1 - p);
        nc += p < 0.5;
      }
      ++nt;
    }
  }
  out.add("loglik", num(ll));
  out.add("prob", mnode(P));
  out.add("correct", ival(nc));
  out.add("total", ival(nt));
  out.add("accuracy", num(nt > 0 ? static_cast<double>(nc) / nt : 0.0));
  return true;
}

bool k_irtpst(const Node &a, Node &out, std::string &err) {
  const Node *c = a.at("chain");
  bool st;
  if (!getb(a, "std", st, err)) return false;
  if (!c || c->k != Node::NUM || c->dim.size() < 2 || c->dim.size() > 3 || c->d.empty()) {
    err = "`chain` must be a numeric array: draws by legislators (by dimensions).";
    return false;
  }
  const int S = c->dim[0], n = c->dim[1], d = c->dim.size() == 3 ? c->dim[2] : 1;
  std::vector<double> ch = c->d;
  for (double e : ch) {
    if (!std::isfinite(e)) {
      err = "`chain` must hold finite numbers.";
      return false;
    }
  }
  auto at = [&](int s, int i, int k) -> double & { return ch[z(s) + z(S) * (z(i) + z(n) * z(k))]; };
  if (st) {
    for (int s = 0; s < S; ++s) {
      for (int k = 0; k < d; ++k) {
        std::vector<double> col(z(n));
        for (int i = 0; i < n; ++i) col[z(i)] = at(s, i, k);
        const double m = mean(col);
        double sd = sdv(col.data(), n);
        if (!(sd > 0)) sd = 1;
        for (int i = 0; i < n; ++i) at(s, i, k) = (at(s, i, k) - m) / sd;
      }
    }
  }
  Mat mn(n, d), sd(n, d), lo(n, d), hi(n, d);
  for (int i = 0; i < n; ++i) {
    for (int k = 0; k < d; ++k) {
      std::vector<double> v(z(S));
      for (int s = 0; s < S; ++s) v[z(s)] = at(s, i, k);
      mn(i, k) = mean(v);
      sd(i, k) = sdv(v.data(), S);
      lo(i, k) = quant(v, 0.025);
      hi(i, k) = quant(v, 0.975);
    }
  }
  out.add("mean", mnode(mn));
  out.add("sd", mnode(sd));
  out.add("lo", mnode(lo));
  out.add("hi", mnode(hi));
  out.add("nsamp", ival(S));
  out.add("std", lgl(st));
  return true;
}

// --- ordered Optimal Classification ------------------------------------------

bool k_ordocl(const Node &a, Node &out, std::string &err) {
  Mat Y;
  int nd, mi;
  if (!getm(a, "Y", Y, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 1)) return false;
  const int n = Y.r, m = Y.c;
  if (!dims(static_cast<double>(n) + m, nd, err)) return false;
  rng g(42);
  const Mat X = rnormm(g, n, nd);
  std::vector<std::vector<double>> cats(z(m)), cuts(z(m));
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(Y(i, j))) cats[z(j)].push_back(Y(i, j));
    }
    std::sort(cats[z(j)].begin(), cats[z(j)].end());
    cats[z(j)].erase(std::unique(cats[z(j)].begin(), cats[z(j)].end()), cats[z(j)].end());
    const int nk = static_cast<int>(cats[z(j)].size()) - 1;
    cuts[z(j)] = nk > 0 ? seqn(-1, 1, nk) : std::vector<double>(1, 0.0);
  }
  Mat nr = rnormm(g, m, nd);
  for (int j = 0; j < m; ++j) {
    double l = 0;
    for (int k = 0; k < nd; ++k) l += nr(j, k) * nr(j, k);
    l = std::sqrt(l) + 1e-12;
    for (int k = 0; k < nd; ++k) nr(j, k) /= l;
  }
  int correct = 0, total = 0, it = 0;
  for (it = 1; it <= mi; ++it) {
    if (intr(err, static_cast<double>(n) * m * nd)) return false;
    const int old = correct;
    correct = total = 0;
    std::vector<std::vector<double>> pj(z(m));
    std::vector<std::vector<int>> act(z(m));
    for (int j = 0; j < m; ++j) {
      const std::vector<double> &cj = cats[z(j)];
      if (cj.size() < 2) continue;
      std::vector<double> sc = cuts[z(j)];
      std::sort(sc.begin(), sc.end());
      for (int i = 0; i < n; ++i) {
        if (std::isnan(Y(i, j))) continue;
        double p = 0;
        for (int k = 0; k < nd; ++k) p += X(i, k) * nr(j, k);
        const int ac = static_cast<int>(std::lower_bound(cj.begin(), cj.end(), Y(i, j)) - cj.begin());
        pj[z(j)].push_back(p);
        act[z(j)].push_back(ac);
        correct += fint(p, sc) == ac;
        ++total;
      }
    }
    for (int j = 0; j < m; ++j) {
      const int K = static_cast<int>(cats[z(j)].size());
      if (K < 2) continue;
      const std::vector<double> &p = pj[z(j)];
      const std::vector<int> &ac = act[z(j)];
      std::vector<std::size_t> o(p.size());
      for (std::size_t t = 0; t < o.size(); ++t) o[t] = t;
      std::stable_sort(o.begin(), o.end(), [&](std::size_t u, std::size_t v) { return p[u] < p[v]; });
      std::vector<double> nc(z(K - 1));
      for (int k = 1; k <= K - 1; ++k) {
        std::vector<double> bv;
        for (std::size_t ii = 0; ii + 1 < o.size(); ++ii) {
          if (ac[o[ii]] <= k - 1 && ac[o[ii + 1]] > k - 1) bv.push_back((p[o[ii]] + p[o[ii + 1]]) / 2);
        }
        nc[z(k - 1)] = !bv.empty() ? median(bv) : k <= static_cast<int>(cuts[z(j)].size()) ? cuts[z(j)][z(k - 1)] : 0.0;
      }
      cuts[z(j)] = nc;
    }
    if (total > 0 && correct == old && it > 1) break;
  }
  if (it > mi) it = mi;
  Node cl = list();
  for (int j = 0; j < m; ++j) {
    cl.nm.emplace_back("");
    cl.l.push_back(vec(cuts[z(j)]));
  }
  cl.nm.clear();
  out.add("ideal", mnode(X));
  out.add("cuts", cl);
  out.add("normals", mnode(nr));
  out.add("correct", num(total > 0 ? static_cast<double>(correct) / total : 0.0));
  out.add("iter", ival(it));
  return true;
}

// --- anchoring vignettes -----------------------------------------------------

bool k_anchv(const Node &a, Node &out, std::string &err) {
  std::vector<double> y;
  Mat V;
  int K;
  if (!getv(a, "y", y, err, true) || !getm(a, "V", V, err, true) || !geti(a, "ncat", K, err, 2)) return false;
  const int n = static_cast<int>(y.size()), nv = V.c;
  if (V.r != n) {
    err = msg("`V` must have one row per rating in `y`, %d.", n);
    return false;
  }
  // the thresholds, n by ncat - 1
  if (!cap((K - 1.0) * std::max(n, 1), kutil::maxcell, "`ncat` times the ratings", err)) return false;
  std::vector<double> vm(z(nv));
  for (int j = 0; j < nv; ++j) {
    double s = 0, c = 0;
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(V(i, j))) {
        s += V(i, j);
        c += 1;
      }
    }
    vm[z(j)] = c > 0 ? s / c : NaN;
  }
  std::vector<double> ord(z(nv));
  {
    std::vector<int> o(z(nv));
    for (int j = 0; j < nv; ++j) o[z(j)] = j;
    std::stable_sort(o.begin(), o.end(), [&](int p, int q) {
      if (std::isnan(vm[z(p)])) return false;
      if (std::isnan(vm[z(q)])) return true;
      return vm[z(p)] < vm[z(q)];
    });
    for (int j = 0; j < nv; ++j) ord[z(j)] = o[z(j)] + 1;
  }
  Mat th(n, K - 1);
  for (int i = 0; i < n; ++i) {
    std::vector<double> sv;
    for (int j = 0; j < nv; ++j) {
      if (!std::isnan(V(i, j))) sv.push_back(V(i, j));
    }
    if (sv.size() >= 2) {
      std::sort(sv.begin(), sv.end());
      const int ns = static_cast<int>(sv.size());
      for (int k = 1; k <= K - 1; ++k) {
        int ix = static_cast<int>((k - 1) * ns / (K - 1.0));
        ix = std::min(ix + 1, ns);
        th(i, k - 1) = sv[z(ix - 1)];
      }
    } else {
      const std::vector<double> s = seqn(1, K, K - 1);
      for (int k = 0; k < K - 1; ++k) th(i, k) = s[z(k)];
    }
  }
  std::vector<double> cs(z(n)), dif(z(K - 1));
  for (int i = 0; i < n; ++i) {
    if (std::isnan(y[z(i)])) {
      cs[z(i)] = NaN;
      continue;
    }
    std::vector<double> t(z(K - 1));
    for (int k = 0; k < K - 1; ++k) t[z(k)] = th(i, k);
    std::sort(t.begin(), t.end());
    cs[z(i)] = fint(y[z(i)], t);
  }
  for (int k = 0; k < K - 1; ++k) dif[z(k)] = sdv(&th.v[z(k) * z(n)], n);
  out.add("scores", vec(cs));
  out.add("thresh", mnode(th));
  out.add("dif", vec(dif));
  out.add("order", ivec(ord));
  out.add("n", ival(n));
  out.add("nvig", ival(nv));
  return true;
}

// --- INDSCAL -----------------------------------------------------------------

bool k_indscl(const Node &a, Node &out, std::string &err) {
  const Node *dl = a.at("D");
  int nd, mi;
  double tol;
  if (!geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 1) || !getd(a, "tol", tol, err, 0, HUGE_VAL)) {
    return false;
  }
  std::vector<Mat> Ds;
  if (dl && dl->k == Node::LIST) {
    for (const Node &e : dl->l) {
      Node t = list();
      t.add("D", e);
      Mat m;
      if (!getm(t, "D", m, err)) return false;
      Ds.push_back(m);
    }
  } else if (dl && dl->k == Node::NUM && dl->dim.size() == 3) {
    const int K = dl->dim[0], r = dl->dim[1], c = dl->dim[2];
    for (int k = 0; k < K; ++k) {
      Mat m(r, c);
      for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) m(i, j) = dl->d[z(k) + z(K) * (z(i) + z(r) * z(j))];
      }
      Ds.push_back(m);
    }
  }
  if (Ds.empty()) {
    err = "`D` must be a list of square dissimilarity matrices.";
    return false;
  }
  const int K = static_cast<int>(Ds.size()), n = Ds[0].r;
  if (!dims(static_cast<double>(n) + K, nd, err)) return false;
  for (const Mat &m : Ds) {
    if (m.r != n || m.c != n || n < 2) {
      err = "`D` must hold square matrices of one size, two rows at least.";
      return false;
    }
    for (double e : m.v) {
      if (!std::isfinite(e)) {
        err = "`D` must hold finite numbers with no missing values.";
        return false;
      }
    }
  }
  if (nd > n) {
    err = msg("`ndim` must be at most the number of stimuli, %d.", n);
    return false;
  }
  Mat Da(n, n);
  for (const Mat &m : Ds) {
    for (std::size_t e = 0; e < Da.v.size(); ++e) Da.v[e] += m.v[e];
  }
  for (double &e : Da.v) e /= K;
  std::vector<double> ev;
  Mat vc;
  eigsym(dcen(Da), ev, vc);
  Mat X(n, nd), W(K, nd, 1.0);
  for (int s = 0; s < nd; ++s) {
    for (int i = 0; i < n; ++i) X(i, s) = vc(i, s) * std::sqrt(std::max(ev[z(s)], 0.0));
  }
  int it = 0;
  for (it = 1; it <= mi; ++it) {
    if (intr(err, static_cast<double>(K) * n * n * nd)) return false;
    const Mat Xo = X;
    for (int k = 0; k < K; ++k) {
      for (int s = 0; s < nd; ++s) {
        double nu = 0, de = 0;
        for (int j = 0; j < n; ++j) {
          for (int i = 0; i < n; ++i) {
            const double df = X(i, s) - X(j, s);
            const double ds = std::sqrt(df * df + 1e-12);
            nu += Ds[z(k)](i, j) * ds;
            de += ds * ds;
          }
        }
        W(k, s) = std::max(nu / (de + 1e-12), 0.01);
      }
    }
    for (int j = 0; j < n; ++j) {
      for (int s = 0; s < nd; ++s) {
        double nu = 0, de = 0;
        for (int k = 0; k < K; ++k) {
          for (int l = 0; l < n; ++l) {
            if (l == j) continue;
            nu += W(k, s) * Ds[z(k)](j, l) * X(l, s);
            de += W(k, s) * W(k, s) + 1e-12;
          }
        }
        if (de > 0 && std::isfinite(nu)) X(j, s) = nu / de;
      }
    }
    double ch = 0, no = 0;
    for (std::size_t e = 0; e < X.v.size(); ++e) {
      ch += (X.v[e] - Xo.v[e]) * (X.v[e] - Xo.v[e]);
      no += Xo.v[e] * Xo.v[e];
    }
    if (std::sqrt(ch) / (std::sqrt(no) + 1e-12) < tol) break;
  }
  if (it > mi) it = mi;
  double ts = 0;
  for (int k = 0; k < K; ++k) {
    Mat Xw = X;
    for (int s = 0; s < nd; ++s) {
      for (int i = 0; i < n; ++i) Xw(i, s) *= std::sqrt(W(k, s));
    }
    const Mat dm = pdist(Xw);
    for (std::size_t e = 0; e < dm.v.size(); ++e) ts += (Ds[z(k)].v[e] - dm.v[e]) * (Ds[z(k)].v[e] - dm.v[e]);
  }
  out.add("config", mnode(X));
  out.add("weights", mnode(W));
  out.add("stress", num(ts));
  out.add("iter", ival(it));
  out.add("nind", ival(K));
  out.add("nstim", ival(n));
  return true;
}

// --- normal vectors and cutting lines ------------------------------------------

bool k_nrmvec(const Node &a, Node &out, std::string &err) {
  Mat X;
  std::vector<double> y;
  if (!getm(a, "x", X, err) || !getv(a, "y", y, err, true)) return false;
  if (static_cast<int>(y.size()) != X.r) {
    err = msg("`y` must have one value per row of `x`, %d.", X.r);
    return false;
  }
  std::vector<double> yv;
  std::vector<int> rows;
  for (int i = 0; i < X.r; ++i) {
    if (!std::isnan(y[z(i)])) {
      rows.push_back(i);
      yv.push_back(y[z(i)]);
    }
  }
  const int m = static_cast<int>(rows.size()), d = X.c;
  Mat A(m, d + 1);
  for (int t = 0; t < m; ++t) {
    A(t, 0) = 1;
    for (int k = 0; k < d; ++k) A(t, k + 1) = X(rows[z(t)], k);
  }
  std::vector<double> b;
  if (!lsq(A, yv, b)) {
    err = "`x` with a constant is rank deficient on the rows where `y` is known.";
    return false;
  }
  std::vector<double> nv(b.begin() + 1, b.end());
  double l = 0;
  for (double e : nv) l += e * e;
  l = std::sqrt(l);
  if (l > 0) {
    for (double &e : nv) e /= l;
  }
  double sr = 0, st = 0;
  const double my = mean(yv);
  for (int t = 0; t < m; ++t) {
    double f = 0;
    for (int k = 0; k <= d; ++k) f += A(t, k) * b[z(k)];
    sr += (yv[z(t)] - f) * (yv[z(t)] - f);
    st += (yv[z(t)] - my) * (yv[z(t)] - my);
  }
  const double ang = nv.size() >= 2 ? std::atan2(nv[1], nv[0]) : 0.0;
  out.add("normal", vec(nv));
  out.add("degrees", num(ang * 180 / 3.141592653589793238462643383280));
  out.add("radians", num(ang));
  out.add("r2", num(st > 0 ? 1 - sr / st : 0.0));
  out.add("coef", vec(b));
  return true;
}

bool k_cutlns(const Node &a, Node &out, std::string &err) {
  Mat N;
  std::vector<double> cp, xl;
  if (!getm(a, "normals", N, err) || !getv(a, "cuts", cp, err) || !getv(a, "xlim", xl, err)) return false;
  if (N.c != 2 || static_cast<int>(cp.size()) != N.r || xl.size() != 2) {
    err = "`normals` must have two columns, `cuts` one value per row and `xlim` two values.";
    return false;
  }
  const int n = N.r;
  Mat E(n, 4), M(n, 2);
  std::vector<double> an(z(n));
  for (int k = 0; k < n; ++k) {
    const double n1 = N(k, 0), n2 = N(k, 1), c = cp[z(k)];
    if (std::fabs(n2) > 1e-10) {
      const double y1 = (c - n1 * xl[0]) / n2, y2 = (c - n1 * xl[1]) / n2;
      E(k, 0) = xl[0];
      E(k, 1) = y1;
      E(k, 2) = xl[1];
      E(k, 3) = y2;
      M(k, 0) = (xl[0] + xl[1]) / 2;
      M(k, 1) = (y1 + y2) / 2;
    } else if (std::fabs(n1) > 1e-10) {
      const double xc = c / n1;
      E(k, 0) = xc;
      E(k, 1) = -10;
      E(k, 2) = xc;
      E(k, 3) = 10;
      M(k, 0) = xc;
      M(k, 1) = 0;
    }
    an[z(k)] = std::atan2(n2, n1) * 180 / 3.141592653589793238462643383280;
  }
  out.add("ends", mnode(E));
  out.add("angles", vec(an));
  out.add("mids", mnode(M));
  out.add("nlines", ival(n));
  return true;
}

// --- DW-NOMINATE and its bootstrap ---------------------------------------------

bool k_dwnom(const Node &a, Node &out, std::string &err) {
  Mat V;
  int nd, mi;
  if (!getm(a, "votes", V, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "maxit", mi, err, 0)) return false;
  const int n = V.r, m = V.c;
  if (!dims(n + 2.0 * m, nd, err)) return false;
  rng g(42);
  Mat X = rnormm(g, n, nd, 0.5);
  std::vector<double> w(z(nd), 1.0 / nd);
  const double beta = 15.0;
  Mat nv = rnormm(g, m, nd), mid(m, nd);
  for (int j = 0; j < m; ++j) {
    double l = 0;
    for (int k = 0; k < nd; ++k) l += nv(j, k) * nv(j, k);
    l = std::sqrt(l) + 1e-12;
    for (int k = 0; k < nd; ++k) nv(j, k) /= l;
  }
  for (int it = 0; it < mi; ++it) {
    if (intr(err, static_cast<double>(n) * m * nd)) return false;
    for (int j = 0; j < m; ++j) {
      int nvc = 0, ny = 0, nn = 0;
      std::vector<double> sy(z(nd)), sn(z(nd));
      for (int i = 0; i < n; ++i) {
        if (std::isnan(V(i, j))) continue;
        ++nvc;
        if (V(i, j) == 1.0) {
          ++ny;
          for (int k = 0; k < nd; ++k) sy[z(k)] += X(i, k);
        } else {
          ++nn;
          for (int k = 0; k < nd; ++k) sn[z(k)] += X(i, k);
        }
      }
      if (nvc < 2) continue;
      std::vector<double> yc(z(nd)), nc(z(nd)), dir(z(nd));
      double dn = 0;
      for (int k = 0; k < nd; ++k) {
        yc[z(k)] = ny > 0 ? sy[z(k)] / ny : mid(j, k);
        nc[z(k)] = nn > 0 ? sn[z(k)] / nn : mid(j, k);
        dir[z(k)] = yc[z(k)] - nc[z(k)];
        dn += dir[z(k)] * dir[z(k)];
      }
      dn = std::sqrt(dn);
      for (int k = 0; k < nd; ++k) {
        if (dn > 1e-10) nv(j, k) = dir[z(k)] / dn;
        mid(j, k) = 0.5 * (yc[z(k)] + nc[z(k)]);
      }
    }
    for (int i = 0; i < n; ++i) {
      int nvc = 0;
      std::vector<double> acc(z(nd));
      for (int j = 0; j < m; ++j) {
        if (std::isnan(V(i, j))) continue;
        for (int k = 0; k < nd; ++k) acc[z(k)] += V(i, j) == 1.0 ? mid(j, k) + 0.5 * nv(j, k) : mid(j, k) - 0.5 * nv(j, k);
        ++nvc;
      }
      if (nvc >= 2) {
        for (int k = 0; k < nd; ++k) X(i, k) = acc[z(k)] / nvc;
      }
    }
    double mx = 0;
    for (int i = 0; i < n; ++i) {
      double s = 0;
      for (int k = 0; k < nd; ++k) s += X(i, k) * X(i, k);
      mx = std::max(mx, std::sqrt(s));
    }
    if (mx > 1.0) {
      for (double &e : X.v) e /= mx;
    }
  }
  double ll = 0;
  int tot = 0, cor = 0;
  std::vector<double> cp(z(m));
  for (int j = 0; j < m; ++j) {
    for (int k = 0; k < nd; ++k) cp[z(j)] += nv(j, k) * mid(j, k);
    for (int i = 0; i < n; ++i) {
      if (std::isnan(V(i, j))) continue;
      double dy = 0, dn = 0;
      for (int k = 0; k < nd; ++k) {
        const double ey = X(i, k) - (mid(j, k) + 0.5 * nv(j, k)), en = X(i, k) - (mid(j, k) - 0.5 * nv(j, k));
        dy += w[z(k)] * ey * ey;
        dn += w[z(k)] * en * en;
      }
      const double u = beta * (std::exp(-0.5 * dy) - std::exp(-0.5 * dn));
      const double p = std::min(std::max(pnorm(u), 1e-10), 1 - 1e-10);
      const double y = V(i, j);
      ll += y * std::log(p) + (1.0 - y) * std::log(1.0 - p);
      const bool py = p > 0.5;
      cor += (py && y == 1.0) || (!py && y == 0.0);
      ++tot;
    }
  }
  out.add("ideal", mnode(X));
  out.add("weights", vec(w));
  out.add("normals", mnode(nv));
  out.add("mids", mnode(mid));
  out.add("cuts", vec(cp));
  out.add("loglik", num(ll));
  out.add("gmp", num(tot > 0 ? static_cast<double>(cor) / tot : 0.0));
  out.add("ndim", ival(nd));
  return true;
}

bool k_nombt(const Node &a, Node &out, std::string &err) {
  Mat V, X, N;
  std::vector<double> cp;
  int nb, seed;
  if (!getm(a, "votes", V, err, true) || !getm(a, "x", X, err) || !getm(a, "normals", N, err) ||
      !getv(a, "cuts", cp, err) || !geti(a, "nboot", nb, err, 1) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const int n = V.r, m = V.c, d = X.c;
  if (X.r != n || N.r != m || N.c != d || static_cast<int>(cp.size()) != m) {
    err = "`x` must have a row per legislator, `normals` and `cuts` a row (value) per vote, and `normals` the "
          "columns of `x`.";
    return false;
  }
  rng g(seed);
  const double beta = 15.0;
  Mat P(n, m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) {
      double p = 0;
      for (int k = 0; k < d; ++k) p += X(i, k) * N(j, k);
      P(i, j) = 1 / (1 + std::exp(-(beta * (p - cp[z(j)]))));
    }
  }
  if (!cap(static_cast<double>(nb) * n * d, kutil::maxcell, "`nboot` times the ideal points", err)) return false;
  std::vector<double> bp(z(nb) * z(n) * z(d));
  for (int b = 0; b < nb; ++b) {
    if (intr(err, 20.0 * n * m * d)) return false;
    Mat S(n, m);
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < n; ++i) S(i, j) = g.unif() < P(i, j) ? 1.0 : 0.0;
    }
    for (std::size_t e = 0; e < S.v.size(); ++e) {
      if (std::isnan(V.v[e])) S.v[e] = NaN;
    }
    Mat Xb = X;
    for (double &e : Xb.v) e += g.norm() * 0.1;
    for (int k = 0; k < 20; ++k) {
      for (int i = 0; i < n; ++i) {
        int nvc = 0;
        for (int j = 0; j < m; ++j) nvc += !std::isnan(S(i, j));
        if (nvc < 2) continue;
        std::vector<double> gr(z(d));
        for (int j = 0; j < m; ++j) {
          if (std::isnan(S(i, j))) continue;
          double p = 0;
          for (int t = 0; t < d; ++t) p += Xb(i, t) * N(j, t);
          const double r = S(i, j) - 1 / (1 + std::exp(-beta * (p - cp[z(j)])));
          for (int t = 0; t < d; ++t) gr[z(t)] += N(j, t) * r;
        }
        for (int t = 0; t < d; ++t) Xb(i, t) += 0.01 * gr[z(t)];
      }
    }
    for (int t = 0; t < d; ++t) {
      for (int i = 0; i < n; ++i) bp[z(b) + z(nb) * (z(i) + z(n) * z(t))] = Xb(i, t);
    }
  }
  Mat se(n, d), mn(n, d);
  for (int t = 0; t < d; ++t) {
    for (int i = 0; i < n; ++i) {
      const double *p = &bp[z(nb) * (z(i) + z(n) * z(t))];
      se(i, t) = sdv(p, nb);
      double s = 0;
      for (int b = 0; b < nb; ++b) s += p[b];
      mn(i, t) = s / nb;
    }
  }
  out.add("se", mnode(se));
  out.add("mean", mnode(mn));
  out.add("nboot", ival(nb));
  return true;
}

// --- EM binary IRT ---------------------------------------------------------------

namespace {

double emdev(const Mat &o, const Mat &n, bool ab) {
  if (ab) {
    double mx = 0;
    for (std::size_t e = 0; e < o.v.size(); ++e) mx = std::max(mx, std::fabs(o.v[e] - n.v[e]));
    return mx;
  }
  double worst = 0;
  for (int d = 0; d < o.c; ++d) {
    double mo = 0, mn = 0;
    for (int i = 0; i < o.r; ++i) {
      mo += o(i, d);
      mn += n(i, d);
    }
    mo /= o.r;
    mn /= o.r;
    double su = 0, sv = 0, uv = 0;
    for (int i = 0; i < o.r; ++i) {
      const double u = o(i, d) - mo, v = n(i, d) - mn;
      su += u * u;
      sv += v * v;
      uv += u * v;
    }
    worst = std::max(worst, su > 0 && sv > 0 ? 1 - uv / std::sqrt(su * sv) : 1.0);
  }
  return worst;
}

}  // namespace

bool emfit(const Mat &Y, int D, int mi, double tol, const double *xp, const double *bp, const Node *start, bool ab,
           Em &r, std::string &err) {
  const int N = Y.r, J = Y.c;
  if (!dims(static_cast<double>(N) + J, D, err)) return false;
  Mat x(N, D), b(J, D);
  std::vector<double> al(z(J));
  if (!start) {
    Mat Zm(N, J);
    for (int j = 0; j < J; ++j) {
      double s = 0;
      for (int i = 0; i < N; ++i) {
        Zm(i, j) = std::isnan(Y(i, j)) ? 0.0 : Y(i, j) == 1 ? 1.0 : -1.0;
        s += Zm(i, j);
      }
      for (int i = 0; i < N; ++i) Zm(i, j) -= s / N;
    }
    const Mat S = mult(Zm, Zm);
    Mat Vm(N, D);
    for (int i = 0; i < N; ++i) {
      for (int d = 0; d < D; ++d) Vm(i, d) = std::pow(i + 1.0, d + 1.0);
    }
    for (int step = 0; step < 50; ++step) {
      Vm = mul(S, Vm);
      for (int d = 0; d < D; ++d) {
        for (int e = 0; e < d; ++e) {
          double s = 0;
          for (int i = 0; i < N; ++i) s += Vm(i, d) * Vm(i, e);
          for (int i = 0; i < N; ++i) Vm(i, d) -= s * Vm(i, e);
        }
        double nr = 0;
        for (int i = 0; i < N; ++i) nr += Vm(i, d) * Vm(i, d);
        nr = std::sqrt(nr);
        for (int i = 0; i < N; ++i) Vm(i, d) = nr > 0 ? Vm(i, d) / nr : 0.0;
      }
    }
    for (std::size_t e = 0; e < Vm.v.size(); ++e) x.v[e] = Vm.v[e] * std::sqrt(static_cast<double>(N));
  } else {
    if (start->k != Node::LIST || start->l.size() < 3) {
      err = "`start` must be a list of the intercepts, the slopes and the ideal points.";
      return false;
    }
    const Node &s0 = start->l[0], &s1 = start->l[1], &s2 = start->l[2];
    if (s0.d.size() != z(J) || s1.d.size() != z(J) * z(D) || s2.d.size() != z(N) * z(D)) {
      err = msg("`start` must hold %d intercepts, %d by %d slopes and %d by %d ideal points.", J, J, D, N, D);
      return false;
    }
    al = s0.d;
    b.v = s1.d;
    x.v = s2.d;
  }
  int it = 0;
  bool cv = false;
  while (it < mi) {
    ++it;
    if (intr(err, static_cast<double>(N) * J * (D + 1))) return false;
    Mat M(N, J), ys(N, J);
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < N; ++i) {
        double mm = al[z(j)];
        for (int d = 0; d < D; ++d) mm += x(i, d) * b(j, d);
        M(i, j) = mm;
        if (std::isnan(Y(i, j))) {
          ys(i, j) = mm;
        } else if (Y(i, j) == 1) {
          ys(i, j) = mm + dnorm(mm) / pnorm(mm);
        } else {
          ys(i, j) = mm - dnorm(mm) / pnorm(-mm);
        }
      }
    }
    Mat x2(N, D + 1);
    for (int i = 0; i < N; ++i) {
      x2(i, 0) = 1;
      for (int d = 0; d < D; ++d) x2(i, d + 1) = x(i, d);
    }
    Mat A1 = tmul(x2, x2), B;
    for (int d = 0; d <= D; ++d) A1(d, d) += 1 / bp[1];
    if (!inv(A1, B)) {
      err = "the M-step of the roll-call parameters is singular.";
      return false;
    }
    Mat rhs = tmul(x2, ys);
    for (double &e : rhs.v) e += bp[0] / bp[1];
    const Mat ab2 = tr(mul(B, rhs));
    std::vector<double> an(z(J));
    Mat bn(J, D);
    for (int j = 0; j < J; ++j) {
      an[z(j)] = ab2(j, 0);
      for (int d = 0; d < D; ++d) bn(j, d) = ab2(j, d + 1);
    }
    std::vector<double> eba(z(D));
    for (int d = 0; d < D; ++d) {
      for (int j = 0; j < J; ++j) eba[z(d)] += bn(j, d) * an[z(j)];
    }
    Mat A2 = tmul(bn, bn), Ai;
    for (int d = 0; d < D; ++d) A2(d, d) += 1 / xp[1];
    if (!inv(A2, Ai)) {
      err = "the M-step of the ideal points is singular.";
      return false;
    }
    Mat r2 = tr(mul(ys, bn));
    for (int d = 0; d < D; ++d) {
      for (int i = 0; i < N; ++i) r2(d, i) += xp[0] / xp[1] - eba[z(d)];
    }
    const Mat xn = tr(mul(Ai, r2));
    if (it > 1) {
      Mat ao(J, 1), anm(J, 1);
      ao.v = al;
      anm.v = an;
      const double dv = std::max(std::max(emdev(x, xn, ab), emdev(ao, anm, ab)), emdev(b, bn, ab));
      cv = dv < tol;
    }
    al = an;
    b = bn;
    x = xn;
    if (cv) break;
  }
  double ll = 0;
  for (int j = 0; j < J; ++j) {
    for (int i = 0; i < N; ++i) {
      if (std::isnan(Y(i, j))) continue;
      double mm = al[z(j)];
      for (int d = 0; d < D; ++d) mm += x(i, d) * b(j, d);
      ll += pnorm(Y(i, j) == 1 ? mm : -mm, true);
    }
  }
  r = {x, b, al, ll, it, cv};
  return true;
}

bool k_irtem(const Node &a, Node &out, std::string &err) {
  Mat Y;
  int D, mi;
  double tol;
  std::vector<double> xp, bp;
  if (!getm(a, "votes", Y, err, true) || !geti(a, "ndim", D, err, 1) || !geti(a, "maxit", mi, err, 1) ||
      !getd(a, "tol", tol, err, 0, HUGE_VAL) || !getv(a, "xprior", xp, err) || !getv(a, "bprior", bp, err)) {
    return false;
  }
  if (xp.size() != 2 || bp.size() != 2 || !(xp[1] > 0) || !(bp[1] > 0)) {
    err = "`xprior` and `bprior` must each be a mean and a positive variance.";
    return false;
  }
  const Node *cv = a.at("conv");
  if (!cv || cv->k != Node::STR || cv->s.size() != 1 || (cv->s[0] != "cor" && cv->s[0] != "abs")) {
    err = "`conv` must be \"cor\" or \"abs\".";
    return false;
  }
  Em r;
  if (!emfit(Y, D, mi, tol, xp.data(), bp.data(), a.at("start"), cv->s[0] == "abs", r, err)) return false;
  out.add("ideal", mnode(r.x));
  out.add("disc", mnode(r.b));
  out.add("diff", vec(r.a));
  out.add("loglik", num(r.ll));
  out.add("iter", ival(r.it));
  out.add("conv", lgl(r.conv));
  return true;
}

// --- bootstrap of AM and Basic Space positions --------------------------------------

bool k_scboot(const Node &a, Node &out, std::string &err) {
  Mat Z;
  int nb, seed;
  if (!getm(a, "Z", Z, err, true) || !geti(a, "nboot", nb, err, 1) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const Node *mt = a.at("method");
  if (!mt || mt->k != Node::STR || mt->s.size() != 1 ||
      (mt->s[0] != "am" && mt->s[0] != "blackbox" && mt->s[0] != "blackbox_t")) {
    err = "`method` must be \"am\", \"blackbox\" or \"blackbox_t\".";
    return false;
  }
  const std::string &me = mt->s[0];
  const int n = Z.r;
  if (!cap(static_cast<double>(nb) * (n + Z.c), kutil::maxcell, "`nboot` times the positions", err)) return false;
  rng g(seed);
  std::vector<std::vector<double>> bs;
  for (int b = 0; b < nb; ++b) {
    if (intr(err, static_cast<double>(n) * Z.c * Z.c)) return false;
    Mat Zb(n, Z.c);
    for (int i = 0; i < n; ++i) {
      const int r = static_cast<int>(g.index(n));
      for (int j = 0; j < Z.c; ++j) Zb(i, j) = Z(r, j);
    }
    Node o = list();
    if (me == "am") {
      if (!am(Zb, o, err)) return false;
      bs.push_back(o.at("zhat")->d);
    } else {
      if (!bb(me == "blackbox" ? Zb : tr(Zb), 2, 8, o, err)) return false;
      const Node *w = o.at("weights");
      bs.push_back(std::vector<double>(w->d.begin(), w->d.begin() + w->dim[0]));
    }
  }
  const std::size_t p = bs[0].size();
  std::vector<double> se(p), mn(p), lo(p), hi(p);
  for (std::size_t j = 0; j < p; ++j) {
    std::vector<double> v(z(nb));
    for (int b = 0; b < nb; ++b) v[z(b)] = bs[z(b)][j];
    se[j] = sdv(v.data(), nb);
    mn[j] = mean(v);
    lo[j] = quant(v, 0.025);
    hi[j] = quant(v, 0.975);
  }
  out.add("se", vec(se));
  out.add("mean", vec(mn));
  out.add("lo", vec(lo));
  out.add("hi", vec(hi));
  out.add("nboot", ival(nb));
  return true;
}

// --- Wordfish ------------------------------------------------------------------------

bool k_wrdfsh(const Node &a, Node &out, std::string &err) {
  Mat T;
  int mi;
  double tol;
  if (!getm(a, "dtm", T, err) || !geti(a, "maxit", mi, err, 1) || !getd(a, "tol", tol, err, 0, HUGE_VAL)) return false;
  const int n = T.r, m = T.c;
  for (double e : T.v) {
    if (e < 0) {
      err = "`dtm` must hold counts from 0.";
      return false;
    }
  }
  if (n < 2) {
    err = "`dtm` must have two documents at least.";
    return false;
  }
  rng g(42);
  std::vector<double> om(z(n)), psi(z(n)), al(z(m)), be(z(m));
  for (double &e : om) e = g.norm() * 0.5;
  double tot = 0;
  for (int i = 0; i < n; ++i) {
    double s = 0;
    for (int j = 0; j < m; ++j) s += T(i, j);
    psi[z(i)] = std::log(s + 1);
    tot += s;
  }
  for (int j = 0; j < m; ++j) {
    double s = 0;
    for (int i = 0; i < n; ++i) s += T(i, j);
    al[z(j)] = std::log(s / tot + 1e-10);
  }
  for (double &e : be) e = g.norm() * 0.1;
  auto clamp = [](double e) { return std::min(std::max(e, -20.0), 20.0); };
  int it = 0;
  for (it = 1; it <= mi; ++it) {
    if (intr(err, static_cast<double>(n) * m)) return false;
    const std::vector<double> oo = om;
    const double mo = mean(om), so = sdv(om.data(), n) + 1e-12;
    std::vector<double> nw(z(n));
    for (int i = 0; i < n; ++i) {
      double gg = 0, h = -1.0;
      for (int j = 0; j < m; ++j) {
        const double mu = std::exp(clamp(psi[z(i)] + al[z(j)] + be[z(j)] * om[z(i)]));
        gg += be[z(j)] * (T(i, j) - mu);
        h -= be[z(j)] * be[z(j)] * mu;
      }
      nw[z(i)] = om[z(i)] - gg / h;
    }
    const double m2 = mean(nw), s2 = sdv(nw.data(), n);
    for (int i = 0; i < n; ++i) om[z(i)] = (nw[z(i)] - m2) / (s2 + 1e-12);
    for (int j = 0; j < m; ++j) {
      al[z(j)] += be[z(j)] * mo;
      be[z(j)] *= so;
    }
    for (int j = 0; j < m; ++j) {
      double ga = 0, ha = -0.01;
      std::vector<double> mu(z(n));
      for (int i = 0; i < n; ++i) {
        mu[z(i)] = std::exp(clamp(psi[z(i)] + al[z(j)] + be[z(j)] * om[z(i)]));
        ga += T(i, j) - mu[z(i)];
        ha -= mu[z(i)];
      }
      al[z(j)] -= ga / ha;
      double gb = 0, hb = -0.01;
      for (int i = 0; i < n; ++i) {
        gb += om[z(i)] * (T(i, j) - mu[z(i)]);
        hb -= om[z(i)] * om[z(i)] * mu[z(i)];
      }
      be[z(j)] -= gb / hb;
    }
    double ch = 0, no = 0;
    for (int i = 0; i < n; ++i) {
      ch += (om[z(i)] - oo[z(i)]) * (om[z(i)] - oo[z(i)]);
      no += oo[z(i)] * oo[z(i)];
    }
    if (std::sqrt(ch) / (std::sqrt(no) + 1e-12) < tol) break;
  }
  if (it > mi) it = mi;
  double ll = 0;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < m; ++j) {
      const double mu = std::exp(clamp(psi[z(i)] + al[z(j)] + be[z(j)] * om[z(i)]));
      ll += T(i, j) * std::log(mu + 1e-15) - mu;
    }
  }
  out.add("pos", vec(om));
  out.add("beta", vec(be));
  out.add("alpha", vec(al));
  out.add("psi", vec(psi));
  out.add("loglik", num(ll));
  out.add("iter", ival(it));
  return true;
}

}  // namespace svk
