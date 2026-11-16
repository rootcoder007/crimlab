// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
// libcrimlab, Bayesian scaling and ideal points: Bayesian Aldrich-McKelvey
// (Hare et al. 2015), Bayesian metric MDS and unfolding (Bakker and Poole
// 2013), the Clinton-Jackman-Rivers binary IRT (Albert and Chib 1993
// Gibbs), alpha-NOMINATE (Carroll et al. 2013), the ordinal factor model
// of Quinn (2004) with the cutpoint steps of Cowles (1996), and the
// dynamic IRT of Martin and Quinn (2002). Conjugate steps where the model
// has them, slice sampling (Neal 2003) elsewhere. The Clinton-Jackman-
// Rivers sweeps draw their uniforms from Philox4x32-10 (Salmon et al.
// 2011), the others from R's generator after set.seed(seed).
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

#include "kutil.h"
#include "svk.h"

namespace svk {

namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();
const double Inf = std::numeric_limits<double>::infinity();

std::size_t z(int i) { return static_cast<std::size_t>(i); }

// What a sampler keeps and how far it counts, checked before anything is
// made: the values kept over all draws within kutil::maxcell, and the
// iterations (burn-in and draws, thinned) within int.
const char *const kept = "`nsamp` times the values kept per draw";
bool sizes(double cells, double iters, const char *what, std::string &err) {
  return cap(cells, kutil::maxcell, kept, err) && cap(iters, 2147483647.0, what, err);
}

bool intr(std::string &err, double w) {
  if (!kutil::stop(w)) return false;
  err.clear();
  return true;
}

double mean(const std::vector<double> &v) {
  double s = 0;
  for (double e : v) s += e;
  return s / static_cast<double>(v.size());
}

// N(mu, sd^2) truncated to [lo, hi] by inversion, one uniform each.
double rtn(rng &g, double mu, double sd, double lo, double hi) {
  const double pl = pnorm((lo - mu) / sd), ph = pnorm((hi - mu) / sd);
  const double u = pl + g.unif() * (ph - pl);
  double o = mu + sd * qnorm(u);
  if (!std::isfinite(o)) o = std::min(std::max(mu, lo), hi);
  return std::min(std::max(o, lo), hi);
}

// The slice sampler of Neal (2003, section 4), stepping out (width w, at
// most m steps) and shrinkage, for k conditionally independent
// coordinates at once; f maps candidates to log densities.
using Lf = std::function<void(const std::vector<double> &, std::vector<double> &)>;

std::vector<double> slice(rng &g, const Lf &f, const std::vector<double> &x0, double w, double lower = -Inf,
                          double upper = Inf, int m = 3) {
  const std::size_t k = x0.size();
  std::vector<double> fy, y(k), L(k), R(k), J(k), K(k), fv;
  f(x0, fy);
  for (std::size_t i = 0; i < k; ++i) y[i] = fy[i] - g.expo();
  for (std::size_t i = 0; i < k; ++i) L[i] = x0[i] - w * g.unif();
  for (std::size_t i = 0; i < k; ++i) {
    R[i] = std::min(L[i] + w, upper);
    L[i] = std::max(L[i], lower);
  }
  for (std::size_t i = 0; i < k; ++i) J[i] = std::floor(m * g.unif());
  for (std::size_t i = 0; i < k; ++i) K[i] = (m - 1) - J[i];
  std::vector<char> go(k);
  for (;;) {
    bool any = false;
    for (std::size_t i = 0; i < k; ++i) any |= (go[i] = J[i] > 0 && L[i] > lower);
    if (!any) break;
    f(L, fv);
    any = false;
    for (std::size_t i = 0; i < k; ++i) {
      if (go[i]) any |= (go[i] = y[i] < fv[i]);
    }
    if (!any) break;
    for (std::size_t i = 0; i < k; ++i) {
      if (!go[i]) continue;
      L[i] = std::max(L[i] - w, lower);
      J[i] -= 1;
    }
  }
  for (;;) {
    bool any = false;
    for (std::size_t i = 0; i < k; ++i) any |= (go[i] = K[i] > 0 && R[i] < upper);
    if (!any) break;
    f(R, fv);
    any = false;
    for (std::size_t i = 0; i < k; ++i) {
      if (go[i]) any |= (go[i] = y[i] < fv[i]);
    }
    if (!any) break;
    for (std::size_t i = 0; i < k; ++i) {
      if (!go[i]) continue;
      R[i] = std::min(R[i] + w, upper);
      K[i] -= 1;
    }
  }
  std::vector<double> x = x0, cand(k);
  std::vector<char> todo(k, 1);
  for (int it = 0; it < 200; ++it) {
    cand = x0;
    for (std::size_t i = 0; i < k; ++i) {
      if (todo[i]) cand[i] = L[i] + g.unif() * (R[i] - L[i]);
    }
    f(cand, fv);
    bool any = false;
    for (std::size_t i = 0; i < k; ++i) {
      if (todo[i] && fv[i] > y[i]) {
        x[i] = cand[i];
        todo[i] = 0;
      }
      any |= todo[i] != 0;
    }
    if (!any) break;
    for (std::size_t i = 0; i < k; ++i) {
      if (!todo[i]) continue;
      if (cand[i] < x0[i]) {
        L[i] = cand[i];
      } else {
        R[i] = cand[i];
      }
    }
  }
  return x;
}

std::vector<double> stdz(const std::vector<double> &v) {
  const double m = mean(v), s = sdv(v.data(), static_cast<int>(v.size()));
  std::vector<double> o(v.size());
  for (std::size_t i = 0; i < v.size(); ++i) o[i] = (v[i] - m) / s;
  return o;
}

// Squared distances between the rows of A and B, dimension by dimension.
Mat sqd(const Mat &A, const Mat &B) {
  Mat o(A.r, B.r);
  for (int k = 0; k < A.c; ++k) {
    for (int j = 0; j < B.r; ++j) {
      for (int i = 0; i < A.r; ++i) {
        const double e = A(i, k) - B(j, k);
        o(i, j) += e * e;
      }
    }
  }
  return o;
}

Mat prot(const Mat &X, const Mat &T) {
  Mat u, v;
  std::vector<double> d;
  svd(tmul(X, T), u, d, v);
  return mult(u, v);
}

std::vector<double> cmean(const Mat &X) {
  std::vector<double> m(z(X.c));
  for (int k = 0; k < X.c; ++k) {
    double s = 0;
    for (int i = 0; i < X.r; ++i) s += X(i, k);
    m[z(k)] = s / X.r;
  }
  return m;
}

Mat ccen(const Mat &X, const std::vector<double> &m) {
  Mat o = X;
  for (int k = 0; k < X.c; ++k) {
    for (int i = 0; i < X.r; ++i) o(i, k) -= m[z(k)];
  }
  return o;
}

// Draws (S, n, d) as one vector; the s-th n by d matrix out and back.
struct Arr3 {
  int S, n, d;
  std::vector<double> v;
  Arr3(int s, int nn, int dd) : S(s), n(nn), d(dd), v(z(s) * z(nn) * z(dd)) {}
  double &at(int s, int i, int k) { return v[z(s) + z(S) * (z(i) + z(n) * z(k))]; }
  Mat get(int s) {
    Mat m(n, d);
    for (int k = 0; k < d; ++k) {
      for (int i = 0; i < n; ++i) m(i, k) = at(s, i, k);
    }
    return m;
  }
  void put(int s, const Mat &m) {
    for (int k = 0; k < d; ++k) {
      for (int i = 0; i < n; ++i) at(s, i, k) = m(i, k);
    }
  }
  Mat mean() {
    Mat m(n, d);
    for (int k = 0; k < d; ++k) {
      for (int i = 0; i < n; ++i) {
        double t = 0;
        for (int s = 0; s < S; ++s) t += at(s, i, k);
        m(i, k) = t / S;
      }
    }
    return m;
  }
  Mat sd() {
    Mat m(n, d);
    for (int k = 0; k < d; ++k) {
      for (int i = 0; i < n; ++i) m(i, k) = sdv(&at(0, i, k), S);
    }
    return m;
  }
  Node node() { return anode(v, S, n, d); }
};

// The truncated gamma draw of the precision under the U(0, 10) prior of
// Bakker and Poole (2013).
double bptau(rng &g, double n, double sse) {
  const double sh = n / 2 + 1, rt = sse / 2;
  return qgamma(g.unif() * pgamma(10, sh, rt), sh, rt);
}

// Classical scaling (Torgerson-Gower) of the start of the MDS sampler.
Mat cmds(const Mat &D, int k) {
  const int m = D.r;
  Mat D2(m, m);
  for (std::size_t e = 0; e < D.v.size(); ++e) D2.v[e] = D.v[e] * D.v[e];
  std::vector<double> rm(z(m));
  double gm = 0;
  for (int i = 0; i < m; ++i) {
    double s = 0;
    for (int j = 0; j < m; ++j) s += D2(i, j);
    rm[z(i)] = s / m;
  }
  for (double e : D2.v) gm += e;
  gm /= static_cast<double>(m) * m;
  Mat B(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) B(i, j) = -0.5 * (D2(i, j) - rm[z(i)] - rm[z(j)] + gm);
  }
  Mat Bs(m, m);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) Bs(i, j) = (B(i, j) + B(j, i)) / 2;
  }
  std::vector<double> ev;
  Mat vc;
  eigsym(Bs, ev, vc);
  Mat X(m, k);
  for (int c = 0; c < k; ++c) {
    const double l = std::sqrt(std::max(ev[z(c)], 0.0));
    for (int i = 0; i < m; ++i) X(i, c) = vc(i, c) * l;
  }
  return X;
}

}  // namespace

// --- Bayesian Aldrich-McKelvey ---------------------------------------------------

bool k_bamscl(const Node &a, Node &out, std::string &err) {
  Mat Z0;
  int S, B, pol, seed;
  if (!getm(a, "Z", Z0, err, true) || !geti(a, "nsamp", S, err, 2) || !geti(a, "burn", B, err, 0) ||
      !geti(a, "polarity", pol, err, 1) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const int q = Z0.c;
  int cols = 0;
  for (int j = 0; j < q; ++j) {
    bool f = false;
    for (int i = 0; i < Z0.r; ++i) f = f || !std::isnan(Z0(i, j));
    cols += f;
  }
  if (cols < 2) {
    err = "Bayesian Aldrich-McKelvey scaling needs placements of at least two stimuli (columns).";
    return false;
  }
  if (pol > q) {
    err = msg("`polarity` must be a column of `Z`, at most %d.", q);
    return false;
  }
  --pol;
  bool pf = false;
  for (int i = 0; i < Z0.r; ++i) pf = pf || !std::isnan(Z0(i, pol));
  if (!pf) {
    err = "The `polarity` stimulus has no placements.";
    return false;
  }
  rng g(seed);
  std::vector<int> rows;
  for (int i = 0; i < Z0.r; ++i) {
    bool f = false;
    for (int j = 0; j < q; ++j) f = f || !std::isnan(Z0(i, j));
    if (f) rows.push_back(i);
  }
  const int N = static_cast<int>(rows.size());
  if (!sizes(2.0 * S * (N + q), static_cast<double>(B) + S, "`burn` plus `nsamp`", err)) return false;
  Mat Z(N, q), ob(N, q);
  for (int t = 0; t < N; ++t) {
    for (int j = 0; j < q; ++j) {
      const double e = Z0(rows[z(t)], j);
      ob(t, j) = std::isnan(e) ? 0.0 : 1.0;
      Z(t, j) = std::isnan(e) ? 0.0 : e;
    }
  }
  std::vector<double> lower(z(q), -100.0), upper(z(q), 100.0);
  upper[z(pol)] = 0;
  std::vector<double> zs(z(q));
  for (int j = 0; j < q; ++j) {
    double s = 0, c = 0;
    for (int i = 0; i < N; ++i) {
      s += Z(i, j) * ob(i, j);
      c += ob(i, j);
    }
    zs[z(j)] = c > 0 ? s / c : NaN;
  }
  {
    std::vector<double> fin;
    for (double e : zs) {
      if (!std::isnan(e)) fin.push_back(e);
    }
    const double m = mean(fin), sd = sdv(fin.data(), static_cast<int>(fin.size()));
    for (double &e : zs) {
      e = (e - m) / sd;
      if (!std::isfinite(e)) e = 0;
    }
  }
  if (zs[z(pol)] > 0) {
    for (double &e : zs) e = -e;
  }
  for (int j = 0; j < q; ++j) zs[z(j)] = std::min(std::max(zs[z(j)], lower[z(j)] + 1e-8), upper[z(j)] - 1e-8);
  std::vector<double> av(z(N)), bv(z(N), 1.0), ti(z(N), 1.0), tj(z(q), 1.0), ni(z(N)), nj(z(q));
  for (int i = 0; i < N; ++i) {
    double s = 0, c = 0;
    for (int j = 0; j < q; ++j) {
      s += Z(i, j) * ob(i, j);
      c += ob(i, j);
    }
    av[z(i)] = s / c;
    ni[z(i)] = c;
  }
  for (int j = 0; j < q; ++j) {
    for (int i = 0; i < N; ++i) nj[z(j)] += ob(i, j);
  }
  double ga = 1, gb = 1;
  Mat oz(S, q), oa(S, N), obb(S, N), ot(S, q);
  int s = 0;
  for (int it = 1; it <= B + S; ++it) {
    if (intr(err, static_cast<double>(N) * q * q * 20)) return false;
    const std::vector<double> zh = stdz(zs);
    for (int i = 0; i < N; ++i) {
      double sw = 0, nu = 0;
      for (int j = 0; j < q; ++j) {
        const double w = ti[z(i)] * tj[z(j)] * ob(i, j);
        sw += w;
        nu += w * (Z(i, j) - bv[z(i)] * zh[z(j)]);
      }
      av[z(i)] = rtn(g, nu / sw, 1 / std::sqrt(sw), -100, 100);
    }
    for (int i = 0; i < N; ++i) {
      double sw = 0, nu = 0;
      for (int j = 0; j < q; ++j) {
        const double w = ti[z(i)] * tj[z(j)] * ob(i, j);
        sw += w * zh[z(j)] * zh[z(j)];
        nu += w * zh[z(j)] * (Z(i, j) - av[z(i)]);
      }
      bv[z(i)] = rtn(g, nu / sw, 1 / std::sqrt(sw), -100, 100);
    }
    Mat R2(N, q);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < q; ++j) {
        const double e = Z(i, j) - av[z(i)] - bv[z(i)] * zh[z(j)];
        R2(i, j) = e * e * ob(i, j);
      }
    }
    for (int j = 0; j < q; ++j) {
      double cs = 0;
      for (int i = 0; i < N; ++i) cs += ti[z(i)] * R2(i, j);
      tj[z(j)] = g.gamma(0.1 + nj[z(j)] / 2, 1 / (0.1 + cs / 2));
    }
    for (int i = 0; i < N; ++i) {
      double rs = 0;
      for (int j = 0; j < q; ++j) rs += R2(i, j) * tj[z(j)];
      ti[z(i)] = g.gamma(ga + ni[z(i)] / 2, 1 / (gb + rs / 2));
    }
    double st = 0, slt = 0;
    for (double e : ti) st += e;
    gb = g.gamma(0.1 + N * ga, 1 / (0.1 + st));
    for (double e : ti) slt += std::log(e);
    const double lgb = std::log(gb);
    ga = slice(
        g,
        [&](const std::vector<double> &v, std::vector<double> &o) {
          o.resize(v.size());
          for (std::size_t t = 0; t < v.size(); ++t) {
            const double x = v[t];
            o[t] = -0.9 * std::log(x) - 0.1 * x + N * x * lgb - N * std::lgamma(x) + (x - 1) * slt;
          }
        },
        {ga}, 8, 0)[0];
    Mat lsd(N, q);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < q; ++j) lsd(i, j) = 0.5 * std::log(ti[z(i)] * tj[z(j)]);
    }
    for (int j = 0; j < q; ++j) {
      zs[z(j)] = slice(
          g,
          [&](const std::vector<double> &v, std::vector<double> &o) {
            o.resize(v.size());
            for (std::size_t t = 0; t < v.size(); ++t) {
              std::vector<double> zz = zs;
              zz[z(j)] = v[t];
              const std::vector<double> sz = stdz(zz);
              double l = dnorm(v[t], true);
              for (int jj = 0; jj < q; ++jj) {
                for (int i = 0; i < N; ++i) {
                  if (ob(i, jj) == 0) continue;
                  const double m = av[z(i)] + bv[z(i)] * sz[z(jj)], e = Z(i, jj) - m;
                  l += lsd(i, jj) - 0.5 * std::exp(2 * lsd(i, jj)) * e * e;
                }
              }
              o[t] = l;
            }
          },
          {zs[z(j)]}, 8, lower[z(j)], upper[z(j)])[0];
    }
    if (it > B) {
      const std::vector<double> sz = stdz(zs);
      for (int j = 0; j < q; ++j) {
        oz(s, j) = sz[z(j)];
        ot(s, j) = tj[z(j)];
      }
      for (int i = 0; i < N; ++i) {
        oa(s, i) = av[z(i)];
        obb(s, i) = bv[z(i)];
      }
      ++s;
    }
  }
  std::vector<double> zm(z(q)), zsd(z(q)), am(z(N)), bm(z(N)), tm(z(q));
  Mat iv(2, q);
  for (int j = 0; j < q; ++j) {
    std::vector<double> c(oz.v.begin() + z(j) * z(S), oz.v.begin() + z(j + 1) * z(S));
    zm[z(j)] = mean(c);
    zsd[z(j)] = sdv(c.data(), S);
    iv(0, j) = quant(c, 0.025);
    iv(1, j) = quant(c, 0.975);
    tm[z(j)] = mean(std::vector<double>(ot.v.begin() + z(j) * z(S), ot.v.begin() + z(j + 1) * z(S)));
  }
  for (int i = 0; i < N; ++i) {
    am[z(i)] = mean(std::vector<double>(oa.v.begin() + z(i) * z(S), oa.v.begin() + z(i + 1) * z(S)));
    bm[z(i)] = mean(std::vector<double>(obb.v.begin() + z(i) * z(S), obb.v.begin() + z(i + 1) * z(S)));
  }
  out.add("zeta", vec(zm));
  out.add("sd", vec(zsd));
  out.add("interval", mnode(iv));
  out.add("a", vec(am));
  out.add("b", vec(bm));
  out.add("tau", vec(tm));
  out.add("draws", mnode(oz));
  out.add("nsamp", ival(S));
  return true;
}

// --- Bayesian metric MDS -----------------------------------------------------------

bool k_bmdsf(const Node &a, Node &out, std::string &err) {
  Mat D;
  int nd, S, B, seed;
  double sg;
  if (!getm(a, "D", D, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "nsamp", S, err, 2) ||
      !geti(a, "burn", B, err, 0) || !getd(a, "sigma", sg, err, 1e-300, HUGE_VAL) ||
      !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  if (D.r != D.c) {
    err = "`D` must be a square dissimilarity matrix.";
    return false;
  }
  const int m = D.r;
  if (nd >= m) {
    err = msg("`ndim` must be below the number of objects, %d.", m);
    return false;
  }
  bool pos = false;
  for (int j = 0; j < m; ++j) {
    for (int i = j + 1; i < m; ++i) pos = pos || (std::isfinite(D(i, j)) && D(i, j) > 0);
  }
  if (!pos) {
    err = "`D` must contain positive distances.";
    return false;
  }
  if (!sizes(static_cast<double>(S) * (static_cast<double>(m) * nd + 1.0), static_cast<double>(B) + S, "`burn` plus `nsamp`", err)) {
    return false;
  }
  rng g(seed);
  std::vector<char> ok(z(m) * z(m));
  Mat LD(m, m), D0(m, m);
  double nobs = 0;
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      D0(i, j) = std::isfinite(D(i, j)) ? D(i, j) : 0.0;
      const int lo = std::min(i, j), hi = std::max(i, j);
      ok[z(j) * z(m) + z(i)] = lo != hi && std::isfinite(D(lo, hi)) && D(lo, hi) > 0;
    }
  }
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < m; ++i) {
      if (ok[z(j) * z(m) + z(i)]) {
        LD(i, j) = std::log(D(i, j));
        if (i < j) nobs += 1;
      }
    }
  }
  Mat X = cmds(D0, nd);
  for (double &e : X.v) e += 0 + 1e-3 * g.norm();
  double tau = 1 / (sg * sg);
  Arr3 keep(S, m, nd);
  std::vector<double> kt(z(S));
  for (int it = 1; it <= B + S; ++it) {
    if (intr(err, static_cast<double>(m) * m * nd * 10)) return false;
    for (int i = 0; i < m; ++i) {
      std::vector<int> oi;
      for (int l = 0; l < m; ++l) {
        if (ok[z(l) * z(m) + z(i)]) oi.push_back(l);
      }
      for (int k = 0; k < nd; ++k) {
        std::vector<double> rest(oi.size());
        for (std::size_t t = 0; t < oi.size(); ++t) {
          double s = 0;
          for (int d = 0; d < nd; ++d) {
            if (d == k) continue;
            const double e = X(oi[t], d) - X(i, d);
            s += e * e;
          }
          rest[t] = s;
        }
        X(i, k) = slice(
            g,
            [&](const std::vector<double> &v, std::vector<double> &o) {
              o.resize(v.size());
              for (std::size_t c = 0; c < v.size(); ++c) {
                double s = 0;
                for (std::size_t t = 0; t < oi.size(); ++t) {
                  const double e = X(oi[t], k) - v[c];
                  const double lh = 0.5 * std::log(std::max(rest[t] + e * e, 1e-24));
                  const double r = LD(i, oi[t]) - lh;
                  s += r * r;
                }
                o[c] = -0.5 * tau * s - v[c] * v[c] / 200;
              }
            },
            {X(i, k)}, 1)[0];
      }
    }
    const Mat sq = sqd(X, X);
    double sse = 0;
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < m; ++i) {
        if (!ok[z(j) * z(m) + z(i)]) continue;
        const double r = LD(i, j) - std::log(std::max(std::sqrt(sq(i, j)), 1e-12));
        sse += r * r;
      }
    }
    tau = bptau(g, nobs, sse / 2);
    if (it > B) {
      keep.put(it - B - 1, X);
      kt[z(it - B - 1)] = tau;
    }
  }
  Mat tg = keep.get(S - 1);
  for (int pass = 0; pass < 2; ++pass) {
    for (int s = 0; s < S; ++s) {
      const Mat Xs = keep.get(s);
      const std::vector<double> mx = cmean(Xs), mt = cmean(tg);
      const Mat Q = prot(ccen(Xs, mx), ccen(tg, mt));
      Mat A = mul(ccen(Xs, mx), Q);
      for (int k = 0; k < nd; ++k) {
        for (int i = 0; i < m; ++i) A(i, k) += mt[z(k)];
      }
      keep.put(s, A);
    }
    tg = keep.mean();
  }
  Mat dm(m, m);
  for (int s = 0; s < S; ++s) {
    const Mat Xs = keep.get(s);
    const Mat sq = sqd(Xs, Xs);
    for (std::size_t e = 0; e < dm.v.size(); ++e) dm.v[e] += std::sqrt(sq.v[e]);
  }
  for (double &e : dm.v) e /= S;
  double sig = 0, tm = 0;
  for (double t : kt) {
    sig += 1 / std::sqrt(t);
    tm += t;
  }
  out.add("pos", mnode(tg));
  out.add("sd", mnode(keep.sd()));
  out.add("dist", mnode(dm));
  out.add("sigma", num(sig / S));
  out.add("tau", num(tm / S));
  out.add("draws", keep.node());
  out.add("nsamp", ival(S));
  return true;
}

// --- Bayesian unfolding ---------------------------------------------------------------

bool k_bunfld(const Node &a, Node &out, std::string &err) {
  Mat P;
  int nd, S, B, seed;
  if (!getm(a, "D", P, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "nsamp", S, err, 2) ||
      !geti(a, "burn", B, err, 0) || !geti(a, "seed", seed, err, -2147483647)) {
    return false;
  }
  const int n = P.r, m = P.c;
  if (nd > std::min(n, m)) {
    err = msg("`ndim` must be at most %d, the smaller dimension of `D`.", std::min(n, m));
    return false;
  }
  if (!sizes(static_cast<double>(S) * ((n + m) * static_cast<double>(nd) + 1), static_cast<double>(B) + S,
             "`burn` plus `nsamp`", err)) {
    return false;
  }
  Mat ok(n, m), LP(n, m), Pz = P;
  double sp = 0, cp = 0;
  for (std::size_t e = 0; e < P.v.size(); ++e) {
    if (std::isfinite(P.v[e]) && P.v[e] > 0) {
      ok.v[e] = 1;
      LP.v[e] = std::log(P.v[e]);
      sp += P.v[e];
      cp += 1;
    }
  }
  if (!(cp > 0)) {
    err = "`D` must contain positive dissimilarities.";
    return false;
  }
  rng g(seed);
  for (std::size_t e = 0; e < P.v.size(); ++e) {
    if (ok.v[e] == 0) Pz.v[e] = sp / cp;
  }
  for (int j = 0; j < m; ++j) {
    double s = 0;
    for (int i = 0; i < n; ++i) s += Pz(i, j);
    s /= n;
    for (int i = 0; i < n; ++i) Pz(i, j) -= s;
  }
  Mat u, v;
  std::vector<double> d;
  svd(Pz, u, d, v);
  Mat X(n, nd), Zs(m, nd);
  for (int k = 0; k < nd; ++k) {
    const double sq = std::sqrt(d[z(k)]);
    for (int i = 0; i < n; ++i) X(i, k) = u(i, k) * sq / std::sqrt(static_cast<double>(n));
    for (int j = 0; j < m; ++j) Zs(j, k) = v(j, k) * sq / std::sqrt(static_cast<double>(m));
  }
  double tau = 1;
  auto ll = [&](int i, int j, double d2) {
    const double r = LP(i, j) - 0.5 * std::log(std::max(d2, 1e-24));
    return r * r * ok(i, j);
  };
  Arr3 kx(S, n, nd), kz(S, m, nd);
  std::vector<double> kt(z(S));
  Mat ds(n, m);
  for (int it = 1; it <= B + S; ++it) {
    if (intr(err, static_cast<double>(n) * m * nd * 20)) return false;
    for (int k = 0; k < nd; ++k) {
      Mat base = sqd(X, Zs);
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
          const double e = X(i, k) - Zs(j, k);
          base(i, j) -= e * e;
        }
      }
      std::vector<double> x0(z(n));
      for (int i = 0; i < n; ++i) x0[z(i)] = X(i, k);
      const std::vector<double> nx = slice(
          g,
          [&](const std::vector<double> &vv, std::vector<double> &o) {
            o.assign(vv.size(), 0.0);
            for (int i = 0; i < n; ++i) {
              double s = 0;
              for (int j = 0; j < m; ++j) {
                const double e = vv[z(i)] - Zs(j, k);
                s += ll(i, j, base(i, j) + e * e);
              }
              o[z(i)] = -0.5 * tau * s - vv[z(i)] * vv[z(i)] / 200;
            }
          },
          x0, 1);
      for (int i = 0; i < n; ++i) X(i, k) = nx[z(i)];
    }
    for (int k = 0; k < nd; ++k) {
      Mat base = sqd(X, Zs);
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
          const double e = X(i, k) - Zs(j, k);
          base(i, j) -= e * e;
        }
      }
      std::vector<double> z0(z(m));
      for (int j = 0; j < m; ++j) z0[z(j)] = Zs(j, k);
      const std::vector<double> nz = slice(
          g,
          [&](const std::vector<double> &vv, std::vector<double> &o) {
            o.assign(vv.size(), 0.0);
            for (int j = 0; j < m; ++j) {
              double s = 0;
              for (int i = 0; i < n; ++i) {
                const double e = X(i, k) - vv[z(j)];
                s += ll(i, j, base(i, j) + e * e);
              }
              o[z(j)] = -0.5 * tau * s - vv[z(j)] * vv[z(j)] / 200;
            }
          },
          z0, 1);
      for (int j = 0; j < m; ++j) Zs(j, k) = nz[z(j)];
    }
    const Mat sq = sqd(X, Zs);
    double tot = 0;
    for (int j = 0; j < m; ++j) {
      for (int i = 0; i < n; ++i) tot += ll(i, j, sq(i, j));
    }
    tau = bptau(g, cp, tot);
    if (it > B) {
      kx.put(it - B - 1, X);
      kz.put(it - B - 1, Zs);
      kt[z(it - B - 1)] = tau;
      for (std::size_t e = 0; e < ds.v.size(); ++e) ds.v[e] += std::sqrt(sq.v[e]);
    }
  }
  Mat tg = kz.get(S - 1);
  for (int pass = 0; pass < 2; ++pass) {
    for (int s = 0; s < S; ++s) {
      const Mat Zd = kz.get(s), Xd = kx.get(s);
      const std::vector<double> mx = cmean(Zd), mt = cmean(tg);
      const Mat Q = prot(ccen(Zd, mx), ccen(tg, mt));
      Mat A = mul(ccen(Xd, mx), Q), C = mul(ccen(Zd, mx), Q);
      for (int k = 0; k < nd; ++k) {
        for (int i = 0; i < n; ++i) A(i, k) += mt[z(k)];
        for (int j = 0; j < m; ++j) C(j, k) += mt[z(k)];
      }
      kx.put(s, A);
      kz.put(s, C);
    }
    tg = kz.mean();
  }
  for (double &e : ds.v) e /= S;
  double sig = 0, tm = 0;
  for (double t : kt) {
    sig += 1 / std::sqrt(t);
    tm += t;
  }
  out.add("stim", mnode(tg));
  out.add("sd", mnode(kz.sd()));
  out.add("ideal", mnode(kx.mean()));
  out.add("dist", mnode(ds));
  out.add("sigma", num(sig / S));
  out.add("tau", num(tm / S));
  out.add("nsamp", ival(S));
  return true;
}

// --- Clinton-Jackman-Rivers ------------------------------------------------------------

namespace {

// n uniforms in (0, 1) of Philox4x32-10 with the key of seed and the
// stream as the third counter word.
std::vector<double> philox(int n, double seed, double stream) {
  const std::uint32_t M0 = 0xD2511F53u, M1 = 0xCD9E8D57u, W0 = 0x9E3779B9u, W1 = 0xBB67AE85u;
  const double t32 = 4294967296.0;
  const std::uint32_t key0 = static_cast<std::uint32_t>(std::fmod(seed, t32)),
                      key1 = static_cast<std::uint32_t>(std::fmod(std::floor(seed / t32), t32));
  const std::uint32_t st = static_cast<std::uint32_t>(std::fmod(stream, t32));
  std::vector<double> o(z(n));
  for (int blk = 0; blk * 4 < n; ++blk) {
    std::uint32_t c[4] = {static_cast<std::uint32_t>(blk), 0u, st, 0u};
    std::uint32_t k0 = key0, k1 = key1;
    for (int r = 0; r < 10; ++r) {
      const std::uint64_t p0 = static_cast<std::uint64_t>(M0) * c[0], p1 = static_cast<std::uint64_t>(M1) * c[2];
      const std::uint32_t n0 = static_cast<std::uint32_t>(p1 >> 32) ^ c[1] ^ k0;
      const std::uint32_t n2 = static_cast<std::uint32_t>(p0 >> 32) ^ c[3] ^ k1;
      c[0] = n0;
      c[1] = static_cast<std::uint32_t>(p1);
      c[2] = n2;
      c[3] = static_cast<std::uint32_t>(p0);
      if (r < 9) {
        k0 += W0;
        k1 += W1;
      }
    }
    for (int w = 0; w < 4 && blk * 4 + w < n; ++w) o[z(blk * 4 + w)] = c[w] / t32 + 0.5 / t32;
  }
  return o;
}

const double AE[8] = {6.65790464350110377720, 5.46378491116411436990, 1.78482653991729133580,
                      2.96560571828504891230e-1, 2.65321895265761230930e-2, 1.24266094738807843860e-3,
                      2.71155556874348757815e-5, 2.01033439929228813265e-7};
const double AF[8] = {1, 5.99832206555887937690e-1, 1.36929880922735805310e-1, 1.48753612908506148525e-2,
                      7.86869131145613259100e-4, 1.84631831751005468180e-5, 1.42151175831644588870e-7,
                      2.04426310338993978564e-15};

double poly(const double *c, double x) {
  double o = c[7];
  for (int i = 6; i >= 0; --i) o = o * x + c[i];
  return o;
}

// Wichura's AS 241 (PPND16).
double as241(double p) {
  static const double A[8] = {3.3871328727963666080, 1.3314166789178437745e2, 1.9715909503065514427e3,
                              1.3731693765509461125e4, 4.5921953931549871457e4, 6.7265770927008700853e4,
                              3.3430575583588128105e4, 2.5090809287301226727e3};
  static const double B[8] = {1, 4.2313330701600911252e1, 6.8718700749205790830e2, 5.3941960214247511077e3,
                              2.1213794301586595867e4, 3.9307895800092710610e4, 2.8729085735721942674e4,
                              5.2264952788528545610e3};
  static const double C[8] = {1.42343711074968357734, 4.63033784615654529590, 5.76949722146069140550,
                              3.64784832476320460504, 1.27045825245236838258, 2.41780725177450611770e-1,
                              2.27238449892691845833e-2, 7.74545014278341407640e-4};
  static const double Dd[8] = {1, 2.05319162663775882187, 1.67638483018380384940, 6.89767334985100004550e-1,
                               1.48103976427480074590e-1, 1.51986665636164571966e-2, 5.47593808499534494600e-4,
                               1.05075007164441684324e-9};
  const double q = p - 0.5;
  if (std::fabs(q) <= 0.425) {
    const double r = 0.180625 - q * q;
    return q * poly(A, r) / poly(B, r);
  }
  double r = q < 0 ? p : 1 - p;
  r = std::sqrt(-std::log(r));
  const double val = r <= 5 ? poly(C, r - 1.6) / poly(Dd, r - 1.6) : poly(AE, r - 5) / poly(AF, r - 5);
  return q < 0 ? -val : val;
}

}  // namespace

bool k_cjrirt(const Node &a, Node &out, std::string &err) {
  Mat Y;
  int D, S, B;
  double bv, seed;
  bool kc;
  if (!getvt(a, "votes", Y, err) || !geti(a, "ndim", D, err, 1) || !geti(a, "nsamp", S, err, 1) ||
      !geti(a, "burn", B, err, 0) || !getd(a, "bvar", bv, err, 1e-300, HUGE_VAL) ||
      !getd(a, "seed", seed, err, 0, 9007199254740992.0) || !getb(a, "chains", kc, err)) {
    return false;
  }
  if (seed != std::floor(seed)) {
    err = "`seed` must be a whole number from 0.";
    return false;
  }
  const int N = Y.r, J = Y.c;
  if (!sizes(static_cast<double>(S) * ((N + J) * static_cast<double>(D) + J), static_cast<double>(B) + S,
             "`burn` plus `nsamp`", err)) {
    return false;
  }
  Mat x;
  if (has(a, "start")) {
    Mat st;
    if (!getm(a, "start", st, err)) return false;
    if (static_cast<int>(st.v.size()) != N * D) {
      err = msg("`start` must hold %d ideal points.", N * D);
      return false;
    }
    x = Mat(N, D);
    x.v = st.v;
  } else {
    Em r;
    const double xp[2] = {0, 1}, bp[2] = {0, 25};
    if (!emfit(Y, D, 500, 1e-6, xp, bp, nullptr, false, r, err)) return false;
    x = r.x;
  }
  std::vector<double> al(z(J));
  Mat be(J, D);
  Arr3 xs(S, N, D), bs(S, J, D);
  Mat as(S, J);
  for (int t = 0; t < B + S; ++t) {
    if (intr(err, static_cast<double>(N) * J * (D + 1))) return false;
    const std::vector<double> u = philox(N * J + J * (D + 1) + N * D, seed, t);
    Mat ys(N, J);
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < J; ++j) {
        double M = -al[z(j)];
        for (int d = 0; d < D; ++d) M += x(i, d) * be(j, d);
        const bool ob = !std::isnan(Y(i, j)), yea = ob && Y(i, j) == 1;
        const double s = yea ? M : -M, v = u[z(i) * z(J) + z(j)];
        double q;
        if (ob && s < -30) {
          const double w = 1 / (s * s);
          const double lp = -0.5 * s * s - std::log(-s) - 0.5 * std::log(2 * 3.141592653589793238462643383280) +
                            std::log(1 - w + 3 * w * w - 15 * w * w * w + 105 * w * w * w * w);
          const double rr = std::sqrt(-(std::log(v) + lp)) - 5;
          q = -(poly(AE, rr) / poly(AF, rr));
        } else {
          q = as241(ob ? v * pnorm(s) : v);
        }
        ys(i, j) = M + (yea ? -q : q);
      }
    }
    Mat xt(N, D + 1);
    for (int i = 0; i < N; ++i) {
      xt(i, 0) = -1;
      for (int d = 0; d < D; ++d) xt(i, d + 1) = x(i, d);
    }
    Mat A1 = tmul(xt, xt), V, L;
    for (int d = 0; d <= D; ++d) A1(d, d) += 1 / bv;
    if (!inv(A1, V) || !chol(V, L)) {
      err = "the conditional covariance of the roll-call parameters is not positive definite.";
      return false;
    }
    const Mat mab = tr(mul(V, tmul(xt, ys)));  // J x (D + 1)
    for (int j = 0; j < J; ++j) {
      std::vector<double> zz(z(D + 1));
      for (int c = 0; c <= D; ++c) zz[z(c)] = as241(u[z(N) * z(J) + z(j) * z(D + 1) + z(c)]);
      for (int c = 0; c <= D; ++c) {
        double s = mab(j, c);
        for (int k = 0; k <= c; ++k) s += zz[z(k)] * L(k, c);
        if (c == 0) {
          al[z(j)] = s;
        } else {
          be(j, c - 1) = s;
        }
      }
    }
    Mat A2 = tmul(be, be), W, K;
    for (int d = 0; d < D; ++d) A2(d, d) += 1;
    if (!inv(A2, W) || !chol(W, K)) {
      err = "the conditional covariance of the ideal points is not positive definite.";
      return false;
    }
    Mat ya = ys;
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < N; ++i) ya(i, j) += al[z(j)];
    }
    const Mat mx = tr(mul(W, tr(mul(ya, be))));
    for (int i = 0; i < N; ++i) {
      for (int c = 0; c < D; ++c) {
        double s = mx(i, c);
        for (int k = 0; k <= c; ++k) s += as241(u[z(N) * z(J) + z(J) * z(D + 1) + z(i) * z(D) + z(k)]) * K(k, c);
        x(i, c) = s;
      }
    }
    if (t >= B) {
      const int s = t - B;
      xs.put(s, x);
      bs.put(s, be);
      for (int j = 0; j < J; ++j) as(s, j) = al[z(j)];
    }
  }
  std::vector<double> am(z(J));
  for (int j = 0; j < J; ++j) {
    double s = 0;
    for (int t = 0; t < S; ++t) s += as(t, j);
    am[z(j)] = s / S;
  }
  out.add("ideal", mnode(xs.mean()));
  out.add("sd", mnode(S > 1 ? xs.sd() : Mat(N, D)));
  out.add("alpha", vec(am));
  out.add("beta", mnode(bs.mean()));
  out.add("nsamp", ival(S));
  if (kc) {
    out.add("xchain", xs.node());
    out.add("achain", mnode(as));
    out.add("bchain", bs.node());
  }
  return true;
}

// --- alpha-NOMINATE -----------------------------------------------------------------------

namespace {

// solve(rWishart(1, v, solve(S))[, , 1]) as R draws it.
bool riwish(rng &g, double v, const Mat &S, Mat &o) {
  const int p = S.r;
  Mat Si, U;
  if (!inv(S, Si) || !chol(Si, U)) return false;
  Mat Zf(p, p);
  for (int j = 0; j < p; ++j) {
    Zf(j, j) = std::sqrt(g.chisq(v - j));
    for (int i = 0; i < j; ++i) Zf(i, j) = g.norm();
  }
  const Mat T = mul(Zf, U);
  return inv(tmul(T, T), o);
}

struct Anom {
  const Mat &V;
  double beta, alpha;
  // log Pr(vote) of legislator i on roll call j at squared distances.
  double ll(int i, int j, double dy, double dn, double b, double al) const {
    const double v = V(i, j);
    if (std::isnan(v)) return 0.0;
    const double w2 = 0.25;
    const double qd = -0.5 * b * w2 * (dy - dn);
    const double nm = b * (std::exp(-0.5 * w2 * dy) - std::exp(-0.5 * w2 * dn));
    const double u = qd + al * (nm - qd);
    return pnorm(v == 1 ? u : -u, true);
  }
};

}  // namespace

bool k_anomin(const Node &a, Node &out, std::string &err) {
  Mat V0;
  int nd, S, B, seed, th;
  double lop, mv;
  bool fix;
  std::vector<double> pv;
  if (!getm(a, "votes", V0, err, true) || !geti(a, "ndim", nd, err, 1) || !geti(a, "nsamp", S, err, 1) ||
      !geti(a, "burn", B, err, 0) || !geti(a, "seed", seed, err, -2147483647) || !geti(a, "thin", th, err, 1) ||
      !getd(a, "lop", lop, err, 0, 0.5) || !getd(a, "minvotes", mv, err, 0, HUGE_VAL) ||
      !getv(a, "polarity", pv, err) || !getb(a, "fixalpha", fix, err)) {
    return false;
  }
  for (double e : V0.v) {
    if (!std::isnan(e) && e != 0 && e != 1) {
      err = "`votes` must hold 1 (yea), 0 (nay) or NA (missing).";
      return false;
    }
  }
  // at most one dimension per legislator, before the polarity is spread
  // over them
  if (!cap(nd, V0.r, "`ndim`", err)) return false;
  // a copy first: assign(n, pv[0]) reads pv[0] after freeing it (libc++)
  if (pv.size() == 1) pv = std::vector<double>(z(nd), pv[0]);
  if (static_cast<int>(pv.size()) != nd) {
    err = msg("`polarity` must give one legislator per dimension, %d.", nd);
    return false;
  }
  std::vector<int> uv, ul;
  for (int j = 0; j < V0.c; ++j) {
    double y = 0, c = 0;
    for (int i = 0; i < V0.r; ++i) {
      if (std::isnan(V0(i, j))) continue;
      c += 1;
      y += V0(i, j) == 1;
    }
    const double mn = std::min(y, c - y) / std::max(c, 1.0);
    if (mn >= lop && c > 0) uv.push_back(j);
  }
  for (int i = 0; i < V0.r; ++i) {
    double c = 0;
    for (int j : uv) c += !std::isnan(V0(i, j));
    if (c >= mv) ul.push_back(i);
  }
  if (uv.size() < 2 || static_cast<int>(ul.size()) < nd + 2) {
    err = "Too few roll calls or legislators survive the `lop` and `minvotes` screens; lower them or supply more "
          "votes.";
    return false;
  }
  std::vector<int> pol(z(nd));
  for (int k = 0; k < nd; ++k) {
    const auto f = std::find(ul.begin(), ul.end(), static_cast<int>(pv[z(k)]) - 1);
    if (pv[z(k)] != std::floor(pv[z(k)]) || f == ul.end()) {
      err = "`polarity` must index legislators kept by the `minvotes` screen.";
      return false;
    }
    pol[z(k)] = static_cast<int>(f - ul.begin());
  }
  const int n = static_cast<int>(ul.size()), m = static_cast<int>(uv.size()), d = nd;
  if (!sizes(static_cast<double>(S) * ((n + 2.0 * m) * d + 2), B + static_cast<double>(S) * th,
             "`burn` plus `nsamp` times `thin`", err)) {
    return false;
  }
  Mat V(n, m);
  for (int j = 0; j < m; ++j) {
    for (int i = 0; i < n; ++i) V(i, j) = V0(ul[z(i)], uv[z(j)]);
  }
  rng g(seed);
  const int nit = B + S * th;
  Mat X(n, d), Yl(m, d), Nl(m, d);
  for (double &e : X.v) e = -1 + 2 * g.unif();
  for (int k = 0; k < d; ++k) X(pol[z(k)], k) = std::fabs(X(pol[z(k)], k));
  for (double &e : Yl.v) e = -1 + 2 * g.unif();
  for (double &e : Nl.v) e = -1 + 2 * g.unif();
  double beta = 10, alpha = fix ? 1 : 0.7;
  Anom an{V, beta, alpha};
  Arr3 dX(S, n, d), dY(S, m, d), dN(S, m, d);
  std::vector<double> db(z(S)), da(z(S));
  auto quad = [&](const Mat &Zm, const Mat &P, int i, int k, double vk) {
    // row i of rowSums((Z %*% P) * Z) with Z[i, k] replaced by vk
    double s = 0;
    for (int c = 0; c < d; ++c) {
      double zp = 0;
      for (int e = 0; e < d; ++e) zp += (e == k ? vk : Zm(i, e)) * P(e, c);
      s += zp * (c == k ? vk : Zm(i, c));
    }
    return s;
  };
  int s = 0;
  for (int it = 1; it <= nit; ++it) {
    if (intr(err, static_cast<double>(n) * m * d * 30)) return false;
    Mat Sx, Sy, Sn;
    if (!riwish(g, n - 1, tmul(X, X), Sx) || !riwish(g, m - 1, tmul(Yl, Yl), Sy) ||
        !riwish(g, m - 1, tmul(Nl, Nl), Sn)) {
      err = "the cross-product of the coordinates is singular.";
      return false;
    }
    Mat DN = sqd(X, Nl);
    for (int k = 0; k < d; ++k) {
      Mat base = sqd(X, Yl);
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
          const double e = X(i, k) - Yl(j, k);
          base(i, j) -= e * e;
        }
      }
      std::vector<double> y0(z(m));
      for (int j = 0; j < m; ++j) y0[z(j)] = Yl(j, k);
      const std::vector<double> ny = slice(
          g,
          [&](const std::vector<double> &v, std::vector<double> &o) {
            o.assign(v.size(), 0.0);
            for (int j = 0; j < m; ++j) {
              double l = 0;
              for (int i = 0; i < n; ++i) {
                const double e = X(i, k) - v[z(j)];
                l += an.ll(i, j, base(i, j) + e * e, DN(i, j), beta, alpha);
              }
              o[z(j)] = l - quad(Yl, Sy, j, k, v[z(j)]) / 2;
            }
          },
          y0, 8);
      for (int j = 0; j < m; ++j) Yl(j, k) = ny[z(j)];
    }
    Mat DY = sqd(X, Yl);
    for (int k = 0; k < d; ++k) {
      Mat base = sqd(X, Nl);
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
          const double e = X(i, k) - Nl(j, k);
          base(i, j) -= e * e;
        }
      }
      std::vector<double> n0(z(m));
      for (int j = 0; j < m; ++j) n0[z(j)] = Nl(j, k);
      const std::vector<double> nn = slice(
          g,
          [&](const std::vector<double> &v, std::vector<double> &o) {
            o.assign(v.size(), 0.0);
            for (int j = 0; j < m; ++j) {
              double l = 0;
              for (int i = 0; i < n; ++i) {
                const double e = X(i, k) - v[z(j)];
                l += an.ll(i, j, DY(i, j), base(i, j) + e * e, beta, alpha);
              }
              o[z(j)] = l - quad(Nl, Sn, j, k, v[z(j)]) / 2;
            }
          },
          n0, 8);
      for (int j = 0; j < m; ++j) Nl(j, k) = nn[z(j)];
    }
    for (int k = 0; k < d; ++k) {
      Mat bY = sqd(X, Yl), bN = sqd(X, Nl);
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
          const double e1 = X(i, k) - Yl(j, k), e2 = X(i, k) - Nl(j, k);
          bY(i, j) -= e1 * e1;
          bN(i, j) -= e2 * e2;
        }
      }
      std::vector<double> x0(z(n));
      for (int i = 0; i < n; ++i) x0[z(i)] = X(i, k);
      const std::vector<double> nx = slice(
          g,
          [&](const std::vector<double> &v, std::vector<double> &o) {
            o.assign(v.size(), 0.0);
            for (int i = 0; i < n; ++i) {
              double l = 0;
              for (int j = 0; j < m; ++j) {
                const double e1 = v[z(i)] - Yl(j, k), e2 = v[z(i)] - Nl(j, k);
                l += an.ll(i, j, bY(i, j) + e1 * e1, bN(i, j) + e2 * e2, beta, alpha);
              }
              o[z(i)] = l - quad(X, Sx, i, k, v[z(i)]) / 2;
            }
          },
          x0, 8);
      for (int i = 0; i < n; ++i) X(i, k) = nx[z(i)];
    }
    DY = sqd(X, Yl);
    DN = sqd(X, Nl);
    auto tot = [&](double b, double al) {
      double l = 0;
      for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) l += an.ll(i, j, DY(i, j), DN(i, j), b, al);
      }
      return l;
    };
    beta = slice(
        g,
        [&](const std::vector<double> &v, std::vector<double> &o) {
          o.resize(v.size());
          for (std::size_t t = 0; t < v.size(); ++t) o[t] = tot(v[t], alpha);
        },
        {beta}, 8, 0)[0];
    if (!fix) {
      alpha = slice(
          g,
          [&](const std::vector<double> &v, std::vector<double> &o) {
            o.resize(v.size());
            for (std::size_t t = 0; t < v.size(); ++t) o[t] = tot(beta, v[t]);
          },
          {alpha}, 8, 0, 1)[0];
    }
    if (it >= B + th && (it - B) % th == 0) {
      dX.put(s, X);
      dY.put(s, Yl);
      dN.put(s, Nl);
      db[z(s)] = beta;
      da[z(s)] = alpha;
      ++s;
    }
  }
  for (int t = 0; t < S; ++t) {
    const Mat Xs = dX.get(t);
    const std::vector<double> mu = cmean(Xs);
    dX.put(t, ccen(Xs, mu));
    dY.put(t, ccen(dY.get(t), mu));
    dN.put(t, ccen(dN.get(t), mu));
  }
  Mat tg = dX.get(S - 1);
  for (int pass = 0; pass < 2; ++pass) {
    for (int t = 0; t < S; ++t) {
      const Mat Q = prot(dX.get(t), tg);
      dX.put(t, mul(dX.get(t), Q));
      dY.put(t, mul(dY.get(t), Q));
      dN.put(t, mul(dN.get(t), Q));
    }
    tg = dX.mean();
  }
  for (int k = 0; k < d; ++k) {
    if (!(tg(pol[z(k)], k) < 0)) continue;
    for (int t = 0; t < S; ++t) {
      for (int i = 0; i < n; ++i) dX.at(t, i, k) = -dX.at(t, i, k);
      for (int j = 0; j < m; ++j) {
        dY.at(t, j, k) = -dY.at(t, j, k);
        dN.at(t, j, k) = -dN.at(t, j, k);
      }
    }
  }
  Mat iv(1, 2);
  iv(0, 0) = quant(da, 0.025);
  iv(0, 1) = quant(da, 0.975);
  std::vector<double> lu, vu;
  for (int i : ul) lu.push_back(i + 1);
  for (int j : uv) vu.push_back(j + 1);
  Node dr = list();
  dr.add("X", dX.node());
  dr.add("Y", dY.node());
  dr.add("N", dN.node());
  dr.add("beta", vec(db));
  dr.add("alpha", vec(da));
  out.add("ideal", mnode(dX.mean()));
  out.add("sd", mnode(dX.sd()));
  out.add("yea", mnode(dY.mean()));
  out.add("nay", mnode(dN.mean()));
  out.add("alpha", num(mean(da)));
  out.add("interval", vec({iv(0, 0), iv(0, 1)}));
  out.add("beta", num(mean(db)));
  out.add("draws", dr);
  out.add("legis", ivec(lu));
  out.add("votes", ivec(vu));
  out.add("ndim", ival(d));
  return true;
}

// --- ordinal IRT (Quinn 2004) -------------------------------------------------------------

namespace {

// prcomp(M, center = TRUE)$x[, 1:k]: the scores on the leading right
// singular vectors of the column-centred matrix.
Mat pcs(const Mat &M0, int k) {
  Mat M = M0;
  for (int j = 0; j < M.c; ++j) {
    double s = 0;
    for (int i = 0; i < M.r; ++i) s += M(i, j);
    s /= M.r;
    for (int i = 0; i < M.r; ++i) M(i, j) -= s;
  }
  Mat u, v;
  std::vector<double> d;
  svd(M, u, d, v);
  Mat Vk(M.c, k);
  for (int c = 0; c < k; ++c) {
    for (int j = 0; j < M.c; ++j) Vk(j, c) = v(j, c);
  }
  return mul(M, Vk);
}

}  // namespace

bool k_ordirt(const Node &a, Node &out, std::string &err) {
  Mat Y;
  int D, S, B, seed;
  double L0;
  if (!getm(a, "Y", Y, err, true) || !geti(a, "ndim", D, err, 1) || !geti(a, "nsamp", S, err, 2) ||
      !geti(a, "burn", B, err, 0) || !geti(a, "seed", seed, err, -2147483647) ||
      !getd(a, "lprec", L0, err, 0, HUGE_VAL)) {
    return false;
  }
  rng g(seed);
  const int n = Y.r, J = Y.c;
  if (D > std::min(n, J)) {
    err = msg("`ndim` must be at most %d, the smaller dimension of `Y`.", std::min(n, J));
    return false;
  }
  if (!sizes(static_cast<double>(S) * (n * static_cast<double>(D) + J * (D + 1.0)), static_cast<double>(B) + S,
             "`burn` plus `nsamp`", err)) {
    return false;
  }
  std::vector<int> nc(z(J));
  Mat ob(n, J);
  for (int j = 0; j < J; ++j) {
    std::vector<double> lv;
    for (int i = 0; i < n; ++i) {
      if (!std::isnan(Y(i, j))) lv.push_back(Y(i, j));
    }
    std::sort(lv.begin(), lv.end());
    lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
    for (int i = 0; i < n; ++i) {
      if (std::isnan(Y(i, j))) continue;
      ob(i, j) = 1;
      Y(i, j) = static_cast<double>(std::lower_bound(lv.begin(), lv.end(), Y(i, j)) - lv.begin()) + 1;
    }
    nc[z(j)] = static_cast<int>(lv.size());
    if (nc[z(j)] < 2) {
      err = "Every item needs at least two observed categories.";
      return false;
    }
  }
  // gam[j][0 .. k]: -Inf, 0, 0.5, 1, ..., Inf (1-based positions 1 .. k + 1)
  std::vector<std::vector<double>> gam(z(J));
  std::vector<double> tune(z(J)), acc(z(J));
  for (int j = 0; j < J; ++j) {
    const int k = nc[z(j)];
    gam[z(j)].push_back(-Inf);
    gam[z(j)].push_back(0);
    for (int c = 1; c <= k - 2; ++c) gam[z(j)].push_back(c * 0.5);
    gam[z(j)].push_back(Inf);
    tune[z(j)] = 0.05 / k;
  }
  Mat M = Y;
  for (int j = 0; j < J; ++j) {
    double s = 0, c = 0;
    for (int i = 0; i < n; ++i) {
      if (ob(i, j) != 0) {
        s += Y(i, j);
        c += 1;
      }
    }
    for (int i = 0; i < n; ++i) {
      if (ob(i, j) == 0) M(i, j) = s / c;
    }
  }
  for (double &e : M.v) e += 0 + 1e-6 * g.norm();
  Mat phi = pcs(M, D);
  for (int k = 0; k < D; ++k) {
    std::vector<double> c(phi.v.begin() + z(k) * z(n), phi.v.begin() + z(k + 1) * z(n));
    const double mm = mean(c), sd = sdv(c.data(), n);
    for (int i = 0; i < n; ++i) phi(i, k) = (phi(i, k) - mm) / sd;
  }
  Mat Lam(J, D + 1), lo(n, J), hi(n, J), ys(n, J);
  auto bounds = [&]() {
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < n; ++i) {
        if (ob(i, j) == 0) continue;
        const int y = static_cast<int>(Y(i, j));
        lo(i, j) = gam[z(j)][z(y - 1)];
        hi(i, j) = gam[z(j)][z(y)];
      }
    }
  };
  bounds();
  Arr3 kp(S, n, D), kl(S, J, D + 1);
  std::vector<std::vector<double>> gs(z(J));
  for (int j = 0; j < J; ++j) gs[z(j)].assign(z(nc[z(j)] - 1), 0.0);
  for (int it = 1; it <= B + S; ++it) {
    if (intr(err, static_cast<double>(n) * J * (D + 1) * 4)) return false;
    Mat mu(n, J);
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < n; ++i) {
        double s = Lam(j, 0);
        for (int d = 0; d < D; ++d) s += phi(i, d) * Lam(j, d + 1);
        mu(i, j) = s;
      }
    }
    for (int j = 0; j < J; ++j) {
      const int k = nc[z(j)];
      if (k <= 2) continue;
      const std::vector<double> gg = gam[z(j)];
      std::vector<double> gp = gg;
      // 1-based positions 3 .. k are 0-based 2 .. k - 1
      for (int c = 2; c <= k - 1; ++c) gp[z(c)] = rtn(g, gg[z(c)], tune[z(j)], gp[z(c - 1)], gg[z(c + 1)]);
      auto ll = [&](const std::vector<double> &q) {
        double s = 0;
        for (int i = 0; i < n; ++i) {
          if (ob(i, j) == 0) continue;
          const int y = static_cast<int>(Y(i, j));
          s += std::log(std::max(pnorm(q[z(y)] - mu(i, j)) - pnorm(q[z(y - 1)] - mu(i, j)), 1e-300));
        }
        return s;
      };
      double c1 = 0, c2 = 0;
      for (int c = 2; c <= k - 1; ++c) {
        c1 += std::log(pnorm((gg[z(c + 1)] - gg[z(c)]) / tune[z(j)]) - pnorm((gp[z(c - 1)] - gg[z(c)]) / tune[z(j)]));
        c2 += std::log(pnorm((gp[z(c + 1)] - gp[z(c)]) / tune[z(j)]) - pnorm((gg[z(c - 1)] - gp[z(c)]) / tune[z(j)]));
      }
      if (std::log(g.unif()) < ll(gp) - ll(gg) + (c1 - c2)) {
        gam[z(j)] = gp;
        acc[z(j)] += 1;
      }
    }
    if (it <= B && it % 50 == 0) {
      for (int j = 0; j < J; ++j) {
        const double rt = acc[z(j)] / 50;
        tune[z(j)] *= rt < 0.2 ? 0.7 : rt > 0.5 ? 1.4 : 1;
        acc[z(j)] = 0;
      }
    }
    if (it == B) std::fill(acc.begin(), acc.end(), 0.0);
    bounds();
    for (int j = 0; j < J; ++j) {
      for (int i = 0; i < n; ++i) {
        if (ob(i, j) != 0) ys(i, j) = rtn(g, mu(i, j), 1, lo(i, j), hi(i, j));
      }
    }
    for (int j = 0; j < J; ++j) {
      std::vector<int> o;
      for (int i = 0; i < n; ++i) {
        if (ob(i, j) != 0) o.push_back(i);
      }
      const int no = static_cast<int>(o.size());
      Mat Xo(no, D + 1), yo(no, 1);
      for (int t = 0; t < no; ++t) {
        Xo(t, 0) = 1;
        for (int d = 0; d < D; ++d) Xo(t, d + 1) = phi(o[z(t)], d);
        yo(t, 0) = ys(o[z(t)], j);
      }
      Mat A = tmul(Xo, Xo), Vv, R;
      for (int d = 0; d <= D; ++d) A(d, d) += L0;
      if (!inv(A, Vv) || !chol(Vv, R)) {
        err = "the conditional covariance of an item's loadings is not positive definite.";
        return false;
      }
      const Mat m1 = mul(Vv, tmul(Xo, yo));
      std::vector<double> zz(z(D + 1));
      for (double &e : zz) e = g.norm();
      for (int c = 0; c <= D; ++c) {
        double s = m1(c, 0);
        for (int k = 0; k <= c; ++k) s += R(k, c) * zz[z(k)];
        Lam(j, c) = s;
      }
    }
    Mat Bm(J, D), Rr(n, J);
    for (int j = 0; j < J; ++j) {
      for (int d = 0; d < D; ++d) Bm(j, d) = Lam(j, d + 1);
      for (int i = 0; i < n; ++i) Rr(i, j) = ys(i, j) - Lam(j, 0);
    }
    std::vector<int> full, part;
    for (int i = 0; i < n; ++i) {
      double c = 0;
      for (int j = 0; j < J; ++j) c += ob(i, j);
      (c == J ? full : part).push_back(i);
    }
    if (!full.empty()) {
      Mat A = tmul(Bm, Bm), W, R;
      for (int d = 0; d < D; ++d) A(d, d) += 1;
      if (!inv(A, W) || !chol(W, R)) {
        err = "the conditional covariance of the ideal points is not positive definite.";
        return false;
      }
      const int nf = static_cast<int>(full.size());
      Mat Rf(nf, J), Zn(nf, D);
      for (int j = 0; j < J; ++j) {
        for (int t = 0; t < nf; ++t) Rf(t, j) = Rr(full[z(t)], j);
      }
      for (double &e : Zn.v) e = g.norm();
      const Mat P1 = mul(mul(Rf, Bm), W), P2 = mul(Zn, R);
      for (int t = 0; t < nf; ++t) {
        for (int d = 0; d < D; ++d) phi(full[z(t)], d) = P1(t, d) + P2(t, d);
      }
    }
    for (int i : part) {
      std::vector<int> o;
      for (int j = 0; j < J; ++j) {
        if (ob(i, j) != 0) o.push_back(j);
      }
      const int no = static_cast<int>(o.size());
      Mat Bo(no, D), ro(no, 1);
      for (int t = 0; t < no; ++t) {
        for (int d = 0; d < D; ++d) Bo(t, d) = Bm(o[z(t)], d);
        ro(t, 0) = Rr(i, o[z(t)]);
      }
      Mat A = tmul(Bo, Bo), W, R;
      for (int d = 0; d < D; ++d) A(d, d) += 1;
      if (!inv(A, W) || !chol(W, R)) {
        err = "the conditional covariance of the ideal points is not positive definite.";
        return false;
      }
      const Mat m1 = mul(W, tmul(Bo, ro));
      std::vector<double> zz(z(D));
      for (double &e : zz) e = g.norm();
      for (int c = 0; c < D; ++c) {
        double s = m1(c, 0);
        for (int k = 0; k <= c; ++k) s += R(k, c) * zz[z(k)];
        phi(i, c) = s;
      }
    }
    if (it > B) {
      kp.put(it - B - 1, phi);
      kl.put(it - B - 1, Lam);
      for (int j = 0; j < J; ++j) {
        for (int c = 1; c < nc[z(j)]; ++c) gs[z(j)][z(c - 1)] += gam[z(j)][z(c)];
      }
    }
  }
  Mat first = kp.get(0);
  for (int s = 0; s < S; ++s) {
    const Mat P = kp.get(s), Lf = kl.get(s);
    Mat L(J, D);
    for (int j = 0; j < J; ++j) {
      for (int d = 0; d < D; ++d) L(j, d) = Lf(j, d + 1);
    }
    Mat Q;
    if (D == 1) {
      Q = Mat(1, 1, L(0, 0) < 0 ? -1.0 : 1.0);
    } else {
      Q = s == 0 ? eye(D) : prot(P, first);
    }
    kp.put(s, mul(P, Q));
    const Mat LQ = mul(L, Q);
    Mat Ln = Lf;
    for (int j = 0; j < J; ++j) {
      for (int d = 0; d < D; ++d) Ln(j, d + 1) = LQ(j, d);
    }
    kl.put(s, Ln);
  }
  const Mat lm = kl.mean();
  Mat disc(J, D);
  std::vector<double> ic(z(J)), ac(z(J));
  for (int j = 0; j < J; ++j) {
    ic[z(j)] = lm(j, 0);
    for (int d = 0; d < D; ++d) disc(j, d) = lm(j, d + 1);
    ac[z(j)] = acc[z(j)] / S;
  }
  Node cl = list();
  for (int j = 0; j < J; ++j) {
    std::vector<double> c = gs[z(j)];
    for (double &e : c) e /= S;
    cl.l.push_back(vec(c));
  }
  out.add("ideal", mnode(kp.mean()));
  out.add("sd", mnode(kp.sd()));
  out.add("disc", mnode(disc));
  out.add("icept", vec(ic));
  out.add("cuts", cl);
  out.add("accept", vec(ac));
  out.add("nsamp", ival(S));
  return true;
}

// --- dynamic IRT (Martin and Quinn 2002) --------------------------------------------------

bool k_dynirt(const Node &a, Node &out, std::string &err) {
  Mat Y;
  int S, B, seed, th;
  double tau2, c0, d0, e0, E0, a0, A0, b0, B0;
  if (!getm(a, "votes", Y, err, true) || !geti(a, "nsamp", S, err, 1) || !geti(a, "burn", B, err, 0) ||
      !geti(a, "seed", seed, err, -2147483647) || !geti(a, "thin", th, err, 1) ||
      !getd(a, "tau2", tau2, err, 0, HUGE_VAL) || !getd(a, "c0", c0, err, -HUGE_VAL, HUGE_VAL) ||
      !getd(a, "d0", d0, err, -HUGE_VAL, HUGE_VAL) || !getd(a, "e0", e0, err, -HUGE_VAL, HUGE_VAL) ||
      !getd(a, "E0", E0, err, 1e-300, HUGE_VAL) || !getd(a, "a0", a0, err, -HUGE_VAL, HUGE_VAL) ||
      !getd(a, "A0", A0, err, 0, HUGE_VAL) || !getd(a, "b0", b0, err, -HUGE_VAL, HUGE_VAL) ||
      !getd(a, "B0", B0, err, 0, HUGE_VAL)) {
    return false;
  }
  const int N = Y.r, K = Y.c;
  if (!cap(B + static_cast<double>(S) * th, 2147483647.0, "`burn` plus `nsamp` times `thin`", err)) return false;
  for (double &e : Y.v) {
    if (!std::isnan(e) && e < 0) e = 0;
    if (!std::isnan(e) && e != 0 && e != 1) {
      err = "`votes` must be 0/1 (or -1/1) with NA for missing.";
      return false;
    }
  }
  const Node *pn = a.at("period");
  std::vector<int> per(z(K));
  std::vector<double> plev;
  Node plv;
  if (pn && pn->k == Node::STR && static_cast<int>(pn->s.size()) == K) {
    std::vector<std::string> lv = pn->s;
    std::sort(lv.begin(), lv.end());
    lv.erase(std::unique(lv.begin(), lv.end()), lv.end());
    for (int k = 0; k < K; ++k) per[z(k)] = static_cast<int>(std::lower_bound(lv.begin(), lv.end(), pn->s[z(k)]) - lv.begin());
    plv.k = Node::STR;
    plv.s = lv;
  } else if (pn && (pn->k == Node::NUM || pn->k == Node::INT) && static_cast<int>(pn->d.size()) == K) {
    plev = pn->d;
    for (double e : plev) {
      if (std::isnan(e)) {
        err = "`period` must not be NA.";
        return false;
      }
    }
    std::sort(plev.begin(), plev.end());
    plev.erase(std::unique(plev.begin(), plev.end()), plev.end());
    for (int k = 0; k < K; ++k) per[z(k)] = static_cast<int>(std::lower_bound(plev.begin(), plev.end(), pn->d[z(k)]) - plev.begin());
    plv = pn->k == Node::INT ? ivec(plev) : vec(plev);
  } else {
    err = msg("`period` must give one period per roll call (ncol(votes) = %d).", K);
    return false;
  }
  int Tn = 0;
  for (int p : per) Tn = std::max(Tn, p + 1);
  rng g(seed);
  Mat Yc = Y;
  for (double &e : Yc.v) {
    if (std::isnan(e)) e = 0.5;
  }
  const Mat pc1 = pcs(Yc, 1);
  std::vector<double> pc(pc1.v);
  const double psd = sdv(pc.data(), N);
  if (psd > 0) {
    pc = stdz(pc);
  } else {
    pc.assign(z(N), 0.0);
  }
  int anc;
  if (has(a, "anchor")) {
    if (!geti(a, "anchor", anc, err, 1)) return false;
    if (anc > N) {
      err = msg("`anchor` must be a row of `votes`, at most %d.", N);
      return false;
    }
    --anc;
  } else {
    anc = static_cast<int>(std::max_element(pc.begin(), pc.end()) - pc.begin());
  }
  Mat theta(N, Tn);
  for (int t = 0; t < Tn; ++t) {
    for (int i = 0; i < N; ++i) theta(i, t) = pc[z(i)];
  }
  std::vector<double> al(z(K), 0.0), be(z(K), 1.0), t2(z(N), tau2);
  const bool est = c0 > 0 && d0 > 0;
  Mat Sth(N, Tn), Sth2(N, Tn);
  std::vector<double> Sa(z(K)), Sb(z(K)), St2(z(N));
  int m = 0;
  Mat lo(N, K), hi(N, K), Zm(N, K);
  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < N; ++i) {
      const double y = Y(i, k);
      lo(i, k) = !std::isnan(y) && y == 1 ? 0.0 : -Inf;
      hi(i, k) = !std::isnan(y) && y == 0 ? 0.0 : Inf;
    }
  }
  for (int it = 1; it <= B + S * th; ++it) {
    if (intr(err, static_cast<double>(N) * K * 4)) return false;
    Mat mu(N, K), pl(N, K), ph(N, K);
    for (int k = 0; k < K; ++k) {
      for (int i = 0; i < N; ++i) {
        mu(i, k) = theta(i, per[z(k)]) * be[z(k)] - al[z(k)];
        pl(i, k) = pnorm(lo(i, k) - mu(i, k));
        ph(i, k) = pnorm(hi(i, k) - mu(i, k));
      }
    }
    for (int k = 0; k < K; ++k) {
      for (int i = 0; i < N; ++i) {
        const double u = g.unif();
        const double p = std::min(std::max(pl(i, k) + u * (ph(i, k) - pl(i, k)), 1e-12), 1 - 1e-12);
        Zm(i, k) = mu(i, k) + qnorm(p);
      }
    }
    for (int k = 0; k < K; ++k) {
      double sx = 0, sx2 = 0, sz = 0, sxz = 0;
      for (int i = 0; i < N; ++i) {
        const double x = theta(i, per[z(k)]);
        sx += x;
        sx2 += x * x;
        sz += Zm(i, k);
        sxz += x * Zm(i, k);
      }
      Mat XtX(2, 2), Vv, R;
      XtX(0, 0) = N + A0;
      XtX(0, 1) = XtX(1, 0) = -sx;
      XtX(1, 1) = sx2 + B0;
      if (!inv(XtX, Vv) || !chol(Vv, R)) {
        err = "the conditional covariance of a roll call's parameters is not positive definite.";
        return false;
      }
      const double z0 = -sz + A0 * a0, z1 = sxz + B0 * b0;
      const double n0 = g.norm(), n1 = g.norm();
      al[z(k)] = Vv(0, 0) * z0 + Vv(0, 1) * z1 + R(0, 0) * n0;
      be[z(k)] = Vv(1, 0) * z0 + Vv(1, 1) * z1 + R(0, 1) * n0 + R(1, 1) * n1;
    }
    Mat po(N, Tn), lin(N, Tn);
    for (int t = 0; t < Tn; ++t) {
      double pr = 0;
      for (int k = 0; k < K; ++k) {
        if (per[z(k)] == t) pr += be[z(k)] * be[z(k)];
      }
      for (int i = 0; i < N; ++i) {
        double s = 0;
        for (int k = 0; k < K; ++k) {
          if (per[z(k)] == t) s += (Zm(i, k) + al[z(k)]) * be[z(k)];
        }
        po(i, t) = pr;
        lin(i, t) = s;
      }
    }
    Mat mf(N, Tn), Pf(N, Tn);
    for (int i = 0; i < N; ++i) {
      double mp = e0, pp = E0;
      for (int t = 0; t < Tn; ++t) {
        const double Pp = pp + t2[z(i)];
        Pf(i, t) = 1 / (1 / Pp + po(i, t));
        mf(i, t) = Pf(i, t) * (mp / Pp + lin(i, t));
        mp = mf(i, t);
        pp = Pf(i, t);
      }
    }
    for (int i = 0; i < N; ++i) theta(i, Tn - 1) = mf(i, Tn - 1) + std::sqrt(Pf(i, Tn - 1)) * g.norm();
    for (int t = Tn - 2; t >= 0; --t) {
      for (int i = 0; i < N; ++i) {
        const double G = Pf(i, t) / (Pf(i, t) + t2[z(i)]);
        const double mb = mf(i, t) + G * (theta(i, t + 1) - mf(i, t)), vb = Pf(i, t) * (1 - G);
        theta(i, t) = mb + std::sqrt(vb) * g.norm();
      }
    }
    if (est && Tn > 1) {
      for (int i = 0; i < N; ++i) {
        double ss = 0;
        for (int t = 1; t < Tn; ++t) ss += (theta(i, t) - theta(i, t - 1)) * (theta(i, t) - theta(i, t - 1));
        t2[z(i)] = 1 / g.gamma((c0 + Tn - 1) / 2, 1 / ((d0 + ss) / 2));
      }
    }
    double ma = 0;
    for (int t = 0; t < Tn; ++t) ma += theta(anc, t);
    if (ma / Tn < 0) {
      for (double &e : theta.v) e = -e;
      for (double &e : be) e = -e;
    }
    if (it >= B + th && (it - B) % th == 0) {
      for (std::size_t e = 0; e < theta.v.size(); ++e) {
        Sth.v[e] += theta.v[e];
        Sth2.v[e] += theta.v[e] * theta.v[e];
      }
      for (int k = 0; k < K; ++k) {
        Sa[z(k)] += al[z(k)];
        Sb[z(k)] += be[z(k)];
      }
      for (int i = 0; i < N; ++i) St2[z(i)] += t2[z(i)];
      ++m;
    }
  }
  Mat tm(N, Tn), ts(N, Tn);
  for (std::size_t e = 0; e < tm.v.size(); ++e) {
    tm.v[e] = Sth.v[e] / m;
    ts.v[e] = std::sqrt(std::max(Sth2.v[e] / m - tm.v[e] * tm.v[e], 0.0));
  }
  for (double &e : Sa) e /= m;
  for (double &e : Sb) e /= m;
  for (double &e : St2) e /= m;
  out.add("theta", mnode(tm));
  out.add("sd", mnode(ts));
  out.add("alpha", vec(Sa));
  out.add("beta", vec(Sb));
  out.add("tau2", vec(St2));
  out.add("periods", plv);
  out.add("nsamp", ival(m));
  out.add("anchor", ival(anc + 1));
  return true;
}

}  // namespace svk
